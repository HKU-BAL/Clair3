use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::Args;

use crate::{
    config::Platform,
    io::{
        compression::read_text_lines, interval::IntervalMap, reference::Reference,
        variant::candidate_positions,
    },
    pipeline::call::{CallConfig, calling_region},
    tensor::{
        full_alignment::create_candidate_details,
        pileup::{PileupOptions, create_pileup},
    },
};

const EXTENDED_WINDOW: i64 = 200;
const REGION_SIZE: i64 = 50_000;
const REFERENCE_REGION_SIZE: i64 = 100_000;
const EXTEND_BP: i64 = 100;
const REFERENCE_ALLELE_GAP: f64 = 0.9;

type Genotype = (i32, i32);

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct UnifyRepresentationArgs {
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long = "var_fn")]
    truth: PathBuf,
    #[arg(long = "ref_fn")]
    reference: PathBuf,
    #[arg(long = "candidate_details_fn")]
    candidate_details: Option<PathBuf>,
    #[arg(long = "output_vcf_fn")]
    output: Option<PathBuf>,
    #[arg(long = "sampleName", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long = "ctgStart")]
    start: Option<i64>,
    #[arg(long = "ctgEnd")]
    end: Option<i64>,
    #[arg(long, default_value_t = 100)]
    max_candidates_distance: i64,
    #[arg(long, default_value_t = 10_000)]
    max_calculate_count: usize,
    #[arg(long, default_value_t = 15)]
    partition_size: usize,
    #[arg(long, default_value_t = 0.15)]
    minimum_allele_gap: f64,
    #[arg(long = "bed_fn")]
    bed: Option<PathBuf>,
    #[arg(long = "vcf_fn")]
    known_vcf: Option<PathBuf>,
    #[arg(long = "extend_bed")]
    extend_bed: Option<PathBuf>,
    #[arg(long, default_value_t = 1000, hide = true)]
    subsample_ratio: usize,
    #[arg(long, default_value_t = 0, hide = true)]
    test_pos: i64,
    #[arg(long)]
    chunk_num: Option<usize>,
    #[arg(long)]
    chunk_id: Option<usize>,
    #[arg(long)]
    min_af: Option<f64>,
    #[arg(long = "bam_fn")]
    bam: Option<PathBuf>,
    #[arg(long, default_value = "pypy3", hide = true)]
    pypy: String,
    #[arg(long = "minCoverage", default_value_t = 2)]
    min_coverage: i64,
    #[arg(long = "minMQ", default_value_t = 5)]
    min_mq: u8,
}

type CandidateRow = (i64, i64, BTreeMap<String, String>);

#[derive(Debug, Clone)]
struct Position {
    start: i64,
    end: i64,
    reference: String,
    alternates: Vec<String>,
    genotype: Genotype,
    confident: bool,
    phased: Option<Genotype>,
    depth: i64,
    alt_list: Vec<(String, i64)>,
    reads: BTreeSet<String>,
    variant_haplotypes: BTreeMap<String, BTreeMap<usize, i64>>,
    haplotype_counts: BTreeMap<usize, i64>,
}

impl Position {
    fn new(start: i64, reference: String, alternates: Vec<String>, genotype: Genotype) -> Self {
        let end = start + reference.len() as i64;
        Self {
            start,
            end,
            reference,
            alternates,
            genotype,
            confident: false,
            phased: None,
            depth: 0,
            alt_list: Vec::new(),
            reads: BTreeSet::new(),
            variant_haplotypes: BTreeMap::new(),
            haplotype_counts: BTreeMap::new(),
        }
    }

    fn set_representation(&mut self, reference: String, alternates: Vec<String>) {
        self.end = self.start + reference.len() as i64;
        self.reference = reference;
        self.alternates = alternates;
    }
}

#[derive(Debug, Default)]
struct ReadInfo {
    haplotype: usize,
    position_alleles: BTreeMap<i64, String>,
    sequences: Vec<(i64, i64, String)>,
}

#[derive(Debug, Clone)]
struct ReferenceWindow {
    start: i64,
    end: i64,
    sequence: Vec<u8>,
}

impl ReferenceWindow {
    fn query(&self, start: i64, end: i64) -> String {
        if end <= start {
            return String::new();
        }
        let start = (start - self.start).max(0) as usize;
        let end = (end - self.start).min(self.sequence.len() as i64).max(0) as usize;
        String::from_utf8_lossy(&self.sequence[start..end]).into_owned()
    }
}

fn remove_common_suffix(reference: &str, alternates: &[String]) -> (String, Vec<String>) {
    if alternates.is_empty() || reference.is_empty() || alternates.iter().any(String::is_empty) {
        return (reference.to_string(), alternates.to_vec());
    }
    let minimum = (reference.len() - 1).min(
        alternates
            .iter()
            .map(|allele| allele.len() - 1)
            .min()
            .unwrap_or(0),
    );
    let mut common = reference.len();
    for allele in alternates {
        common = common.min(allele.len());
        common = reference
            .as_bytes()
            .iter()
            .rev()
            .zip(allele.as_bytes().iter().rev())
            .take(common)
            .take_while(|(left, right)| left == right)
            .count();
    }
    if common > minimum {
        return (reference.to_string(), alternates.to_vec());
    }
    (
        reference[..reference.len() - common].to_string(),
        alternates
            .iter()
            .map(|allele| allele[..allele.len() - common].to_string())
            .collect(),
    )
}

fn decode_alt_info(
    reads: &BTreeMap<String, String>,
    reference_base: &str,
    depth: i64,
    minimum_allele_gap: f64,
) -> (String, Vec<String>, Vec<(String, i64)>) {
    let mut counts: Vec<(String, i64)> = reads
        .iter()
        .filter(|(allele, _)| allele.as_str() != "*")
        .map(|(allele, reads)| {
            (
                allele.clone(),
                reads.split(' ').filter(|name| !name.is_empty()).count() as i64,
            )
        })
        .collect();
    counts.sort_by(|left, right| right.1.cmp(&left.1).then(left.0.cmp(&right.0)));
    let mut alternates = Vec::new();
    let mut deletions = Vec::new();
    let mut longest_deletion = "";
    for (allele, count) in counts.iter().take(2) {
        if allele.contains('*') || allele.contains('#') || depth == 0 {
            continue;
        }
        if *count as f64 / depth as f64 >= minimum_allele_gap {
            match allele.as_bytes()[0] {
                b'X' => alternates.push(allele[1..].to_string()),
                b'I' => alternates.push(allele[1..].to_string()),
                b'D' => {
                    if allele.len() - 1 > longest_deletion.len() {
                        longest_deletion = &allele[1..];
                    }
                    deletions.push(format!("{}{}", reference_base, &allele[1..]));
                }
                _ => {}
            }
        }
    }
    let mut reference = reference_base.to_string();
    if !longest_deletion.is_empty() {
        reference.push_str(longest_deletion);
        alternates = alternates
            .into_iter()
            .map(|allele| format!("{allele}{longest_deletion}"))
            .collect();
        let mut normalized_deletions = Vec::new();
        for deletion in deletions {
            let length = deletion.len() - 1;
            normalized_deletions.push(format!(
                "{}{}",
                reference_base,
                &longest_deletion[length.min(longest_deletion.len())..]
            ));
        }
        alternates.extend(normalized_deletions);
    }
    (reference, alternates, counts)
}

fn lock_variant(candidate: &Position, truth: Option<&Position>) -> bool {
    let Some(truth) = truth else {
        return false;
    };
    if candidate.alternates.len() != truth.alternates.len() {
        return false;
    }
    let candidate_alleles: BTreeSet<(String, String)> = candidate
        .alternates
        .iter()
        .map(|alternate| {
            let (reference, mut alternates) =
                remove_common_suffix(&candidate.reference, std::slice::from_ref(alternate));
            (alternates.remove(0), reference)
        })
        .collect();
    truth.alternates.iter().all(|alternate| {
        let (reference, mut alternates) =
            remove_common_suffix(&truth.reference, std::slice::from_ref(alternate));
        candidate_alleles.contains(&(alternates.remove(0), reference))
    })
}

fn phased_genotype(
    candidate: &Position,
    truth: Option<&Position>,
    reference_base: u8,
    minimum_allele_gap: f64,
    platform: Platform,
) -> Option<Genotype> {
    if candidate.variant_haplotypes.is_empty() {
        return None;
    }
    let total: i64 = candidate.haplotype_counts.values().sum();
    if total == 0 {
        return None;
    }
    let mut phased = (-1, -1);
    for (variant, counts) in &candidate.variant_haplotypes {
        if variant == "*" || variant == "#" {
            continue;
        }
        let hap0 = counts.get(&0).copied().unwrap_or(0);
        if platform == Platform::Ilmn {
            let fraction = hap0 as f64 / total as f64;
            if variant != "R"
                && fraction > REFERENCE_ALLELE_GAP
                && truth.is_some_and(|truth| truth.confident && truth.genotype == (1, 1))
            {
                return Some((1, 1));
            }
            continue;
        }
        let hap1 = counts.get(&1).copied().unwrap_or(0);
        let hap2 = counts.get(&2).copied().unwrap_or(0);
        let hap0 = if hap0 > 3 { hap0 } else { 0 };
        let hap1 = if hap1 > 3 { hap1 } else { 0 };
        let hap2 = if hap2 > 3 { hap2 } else { 0 };
        if variant != "R"
            && variant != "#"
            && (hap0 + hap1 + hap2) as f64 / total as f64 > 1.0 - minimum_allele_gap / 2.0
            && truth.is_some_and(|truth| truth.confident && truth.genotype == (1, 1))
        {
            return Some((1, 1));
        }
        let hap1_fraction = candidate
            .haplotype_counts
            .get(&1)
            .filter(|count| **count != 0)
            .map_or(0.0, |count| hap1 as f64 / *count as f64);
        let hap2_fraction = candidate
            .haplotype_counts
            .get(&2)
            .filter(|count| **count != 0)
            .map_or(0.0, |count| hap2 as f64 / *count as f64);
        if variant == "R" {
            if hap1_fraction >= 1.0 - minimum_allele_gap {
                phased.0 = 0;
            }
            if hap2_fraction >= 1.0 - minimum_allele_gap {
                phased.1 = 0;
            }
            continue;
        }
        let (reference, alternate) = match variant.as_bytes()[0] {
            b'X' => (reference_base.to_string(), variant[1..].to_string()),
            b'I' => (reference_base.to_string(), variant[1..].to_string()),
            b'D' => (
                format!("{}{}", reference_base as char, &variant[1..]),
                reference_base.to_string(),
            ),
            _ => continue,
        };
        let Some(index) = candidate.alternates.iter().position(|allele| {
            let (candidate_reference, mut candidate_alternate) =
                remove_common_suffix(&candidate.reference, std::slice::from_ref(allele));
            candidate_reference == reference && candidate_alternate.remove(0) == alternate
        }) else {
            continue;
        };
        if hap1_fraction >= 1.0 - minimum_allele_gap * 2.0 {
            phased.0 = index as i32 + 1;
        }
        if hap2_fraction >= 1.0 - minimum_allele_gap * 2.0 {
            phased.1 = index as i32 + 1;
        }
    }
    (phased.0 >= 0 && phased.1 >= 0).then_some(phased)
}

fn match_alt_base(candidate: &Position, truth: &Position) -> bool {
    let supported: BTreeSet<&str> = candidate
        .alt_list
        .iter()
        .map(|(allele, _)| allele.as_str())
        .collect();
    truth.alternates.iter().any(|alternate| {
        let (reference, mut alternate) =
            remove_common_suffix(&truth.reference, std::slice::from_ref(alternate));
        let alternate = alternate.remove(0);
        if reference.len() == alternate.len() {
            supported.contains(format!("X{alternate}").as_str())
        } else if reference.len() < alternate.len() {
            supported.contains(format!("I{alternate}").as_str())
        } else {
            supported.contains(format!("D{}", &reference[1..]).as_str())
        }
    })
}

fn unique_genotypes(options: &[BTreeSet<Genotype>]) -> Vec<Vec<Genotype>> {
    fn visit(
        index: usize,
        options: &[BTreeSet<Genotype>],
        current: &mut Vec<Genotype>,
        output: &mut Vec<Vec<Genotype>>,
        seen: &mut BTreeSet<(String, String)>,
    ) {
        if index == options.len() {
            let first = current
                .iter()
                .map(|gt| gt.0.to_string())
                .collect::<String>();
            let second = current
                .iter()
                .map(|gt| gt.1.to_string())
                .collect::<String>();
            let reverse = (second.clone(), first.clone());
            if !seen.contains(&(first.clone(), second.clone())) && !seen.contains(&reverse) {
                seen.insert((first.clone(), second.clone()));
                seen.insert(reverse);
                output.push(current.clone());
            }
            return;
        }
        for genotype in &options[index] {
            current.push(*genotype);
            visit(index + 1, options, current, output, seen);
            current.pop();
        }
    }
    let mut output = Vec::new();
    let mut seen = BTreeSet::new();
    visit(0, options, &mut Vec::new(), &mut output, &mut seen);
    output
}

fn genotype_combinations(
    variants: &[Position],
    partition_truths: &[Position],
    candidates: bool,
    truths: &BTreeMap<i64, Position>,
    maximum: usize,
    simplify: bool,
) -> Vec<BTreeSet<Genotype>> {
    let has_multiple_alternates = candidates
        && partition_truths
            .iter()
            .any(|truth| truth.alternates.len() > 1);
    let mut output = Vec::new();
    for variant in variants {
        let mut genotypes = BTreeSet::new();
        if candidates {
            let truth_phased = truths
                .get(&variant.start)
                .is_some_and(|truth| truth.phased.is_some());
            if let Some(phased) = variant.phased {
                if truth_phased || simplify {
                    genotypes.insert(phased);
                } else {
                    genotypes.insert((0, 0));
                    genotypes.insert(phased);
                }
            } else {
                for first in 0..=variant.alternates.len() as i32 {
                    for second in 0..=variant.alternates.len() as i32 {
                        if simplify
                            && !has_multiple_alternates
                            && first != 0
                            && second != 0
                            && first != second
                        {
                            continue;
                        }
                        genotypes.insert((first, second));
                    }
                }
            }
        } else if variant.confident && variant.phased.is_some() {
            genotypes.insert(variant.phased.unwrap());
        } else {
            genotypes.insert((0, 0));
            genotypes.insert(variant.genotype);
            genotypes.insert((variant.genotype.1, variant.genotype.0));
        }
        output.push(genotypes);
    }
    let count = output.iter().fold(1_usize, |count, options| {
        count.saturating_mul(options.len())
    });
    if count > maximum {
        if simplify {
            return vec![BTreeSet::from([(0, 0)]); variants.len()];
        }
        return genotype_combinations(
            variants,
            partition_truths,
            candidates,
            truths,
            maximum,
            true,
        );
    }
    output
}

fn build_sequences(
    variants: &[Position],
    genotypes: &[Genotype],
    reference: &ReferenceWindow,
    mut previous: [i64; 2],
) -> Option<Vec<String>> {
    let mut sequences = [String::new(), String::new()];
    for (variant, genotype) in variants.iter().zip(genotypes) {
        for haplotype in 0..2 {
            let allele_index = if haplotype == 0 {
                genotype.0
            } else {
                genotype.1
            };
            if variant.start < previous[haplotype] {
                if variant.start == previous[haplotype] - 1 && allele_index != 0 {
                    let alternate = variant.alternates.get(allele_index as usize - 1)?;
                    if alternate.len() > variant.reference.len() {
                        sequences[haplotype].push_str(&alternate[1..]);
                        continue;
                    }
                }
                if allele_index != 0 {
                    return None;
                }
                continue;
            }
            sequences[haplotype].push_str(&reference.query(previous[haplotype], variant.start));
            if allele_index == 0 {
                sequences[haplotype].push(variant.reference.as_bytes()[0] as char);
                previous[haplotype] = variant.start + 1;
            } else {
                let alternate = variant.alternates.get(allele_index as usize - 1)?;
                let (trimmed_reference, trimmed_alternate) =
                    remove_common_suffix(&variant.reference, std::slice::from_ref(alternate));
                previous[haplotype] = variant.start + trimmed_reference.len() as i64;
                sequences[haplotype].push_str(&trimmed_alternate[0]);
            }
        }
    }
    for haplotype in 0..2 {
        if previous[haplotype] != reference.end && previous[haplotype] != 0 {
            sequences[haplotype].push_str(&reference.query(previous[haplotype], reference.end));
        }
    }
    Some(vec![sequences[0].clone(), sequences[1].clone()])
}

fn sequence_combinations(
    variants: &[Position],
    combinations: Vec<Vec<Genotype>>,
    reference: &ReferenceWindow,
) -> BTreeMap<Vec<String>, Vec<Vec<Genotype>>> {
    let mut output = BTreeMap::<Vec<String>, Vec<Vec<Genotype>>>::new();
    for genotypes in combinations {
        let Some(mut sequences) = build_sequences(
            variants,
            &genotypes,
            reference,
            [reference.start, reference.start],
        ) else {
            continue;
        };
        sequences.sort();
        sequences.dedup();
        output.entry(sequences).or_default().push(genotypes);
    }
    output
}

fn read_sequences(
    candidates: &[Position],
    reads: &BTreeMap<String, ReadInfo>,
    reference: &ReferenceWindow,
) -> BTreeMap<String, i64> {
    let names: BTreeSet<&String> = candidates
        .iter()
        .flat_map(|candidate| candidate.reads.iter())
        .collect();
    let mut counts = BTreeMap::new();
    for name in names {
        let Some(read) = reads.get(name) else {
            continue;
        };
        let mut sequence = String::new();
        let mut previous_end = 0_usize;
        for (start, end, alternate) in &read.sequences {
            if *end < reference.start || *start > reference.end {
                continue;
            }
            let start_offset = start - reference.start;
            let end_offset = end - reference.start;
            if start_offset >= previous_end as i64 && start_offset >= 0 {
                sequence.push_str(&reference.query(reference.start + previous_end as i64, *start));
                sequence.push_str(alternate);
                previous_end = end_offset.max(0) as usize;
            }
        }
        if previous_end < reference.sequence.len() {
            sequence
                .push_str(&reference.query(reference.start + previous_end as i64, reference.end));
        }
        *counts.entry(sequence).or_default() += 1;
    }
    counts
}

#[derive(Debug)]
struct ReadMatch {
    candidates: Vec<Position>,
    candidate_genotypes: Vec<Genotype>,
    truths: Vec<Position>,
    truth_genotypes: Vec<Genotype>,
    order: (i64, usize, usize, usize),
}

impl ReadMatch {
    fn new(
        candidates: Vec<Position>,
        candidate_genotypes: Vec<Genotype>,
        truths: Vec<Position>,
        truth_genotypes: Vec<Genotype>,
        match_reads: i64,
    ) -> Self {
        let missing = truths
            .iter()
            .zip(&truth_genotypes)
            .filter(|(truth, genotype)| {
                genotype.0 + genotype.1 < truth.genotype.0 + truth.genotype.1
            })
            .count();
        let nonvariants = candidate_genotypes
            .iter()
            .filter(|genotype| genotype.0 + genotype.1 == 0)
            .count();
        let matched = candidate_genotypes.len() - nonvariants;
        Self {
            candidates,
            candidate_genotypes,
            truths,
            truth_genotypes,
            order: (match_reads, missing, nonvariants, matched),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum VariantKind {
    Candidate,
    Truth,
}

struct Unification<'a> {
    contig: &'a str,
    reference: &'a [u8],
    reference_start: i64,
    partition_size: usize,
    max_candidates_distance: i64,
    max_calculate_count: usize,
}

impl Unification<'_> {
    fn reference_window(
        &self,
        candidates: &[Position],
        truths: &[Position],
    ) -> Result<ReferenceWindow> {
        let start = candidates
            .iter()
            .chain(truths)
            .map(|variant| variant.start)
            .min()
            .context("empty variant partition")?;
        let end = candidates
            .iter()
            .chain(truths)
            .map(|variant| variant.end)
            .max()
            .context("empty variant partition")?;
        let sequence_start = start - self.reference_start - 1;
        let sequence_end = (end + 50 - self.reference_start).min(self.reference.len() as i64);
        if sequence_start < 0 || sequence_start >= sequence_end {
            bail!("variant partition exceeds the loaded reference region");
        }
        Ok(ReferenceWindow {
            start: start - 1,
            end: start - 1 + (sequence_end - sequence_start),
            sequence: self.reference[sequence_start as usize..sequence_end as usize].to_vec(),
        })
    }

    fn find_match_pairs(
        &self,
        candidates: &[Position],
        truths: &[Position],
        reference: &ReferenceWindow,
        truth_positions: &BTreeMap<i64, Position>,
        reads: &BTreeMap<String, ReadInfo>,
    ) -> Option<ReadMatch> {
        if candidates.is_empty() || truths.is_empty() {
            return None;
        }
        let confident = candidates.len() == truths.len()
            && truths.iter().all(|truth| {
                candidates.iter().any(|candidate| {
                    candidate.start == truth.start
                        && candidate.reference == truth.reference
                        && candidate.alternates == truth.alternates
                })
            });
        if confident {
            let genotypes: Vec<Genotype> = truths.iter().map(|truth| truth.genotype).collect();
            return Some(ReadMatch::new(
                truths.to_vec(),
                genotypes.clone(),
                truths.to_vec(),
                genotypes,
                100,
            ));
        }
        let truth_combinations = unique_genotypes(&genotype_combinations(
            truths,
            truths,
            false,
            truth_positions,
            self.max_calculate_count,
            false,
        ));
        let candidate_combinations = unique_genotypes(&genotype_combinations(
            candidates,
            truths,
            true,
            truth_positions,
            self.max_calculate_count,
            false,
        ));
        if truth_combinations
            .len()
            .saturating_mul(candidate_combinations.len())
            > self.max_calculate_count
        {
            return None;
        }
        let truth_sequences = sequence_combinations(truths, truth_combinations, reference);
        let candidate_sequences =
            sequence_combinations(candidates, candidate_combinations, reference);
        let read_sequences = read_sequences(candidates, reads, reference);
        let mut matches = Vec::new();
        for (sequence, candidate_genotypes) in candidate_sequences {
            let Some(truth_genotypes) = truth_sequences.get(&sequence) else {
                continue;
            };
            let support = -sequence
                .iter()
                .map(|sequence| read_sequences.get(sequence).copied().unwrap_or(0))
                .sum::<i64>();
            for genotype in candidate_genotypes {
                matches.push(ReadMatch::new(
                    candidates.to_vec(),
                    genotype,
                    truths.to_vec(),
                    truth_genotypes[0].clone(),
                    support,
                ));
            }
        }
        matches.into_iter().min_by_key(|matched| matched.order)
    }

    fn partitions(
        &self,
        candidates: &[Position],
        truths: &[Position],
        alt_positions: &BTreeMap<i64, Position>,
        truth_positions: &BTreeMap<i64, Position>,
    ) -> Vec<(Vec<Position>, Vec<Position>)> {
        let combination_count = |position: i64, kind: VariantKind| -> usize {
            let alternates = match kind {
                VariantKind::Candidate => alt_positions
                    .get(&position)
                    .map_or(0, |variant| variant.alternates.len()),
                VariantKind::Truth => truth_positions
                    .get(&position)
                    .map_or(0, |variant| variant.alternates.len()),
            };
            if alt_positions
                .get(&position)
                .is_some_and(|variant| variant.phased.is_some())
                && (kind == VariantKind::Candidate || truth_positions.contains_key(&position))
            {
                1
            } else if kind == VariantKind::Candidate {
                (alternates + 1) * alternates / 2
            } else {
                alternates
            }
        };
        let mut variants: Vec<(i64, VariantKind)> = candidates
            .iter()
            .map(|variant| (variant.start, VariantKind::Candidate))
            .chain(
                truths
                    .iter()
                    .map(|variant| (variant.start, VariantKind::Truth)),
            )
            .collect();
        variants.sort_by_key(|variant| *variant);
        let truth_positions_set: BTreeSet<i64> = truths.iter().map(|truth| truth.start).collect();
        let get_variant = |position: i64, kind: VariantKind| -> Option<Position> {
            match kind {
                VariantKind::Candidate => alt_positions.get(&position).cloned(),
                VariantKind::Truth => truth_positions.get(&position).cloned(),
            }
        };
        let mut partitions: Vec<Vec<(i64, VariantKind)>> = Vec::new();
        let mut partition: Vec<(i64, VariantKind)> = Vec::new();
        let mut product = 1_usize;
        for (index, variant) in variants.iter().copied().enumerate() {
            let new_count = product.saturating_mul(combination_count(variant.0, variant.1));
            let matches = if partition.is_empty() {
                true
            } else {
                let type_count = partition.iter().filter(|item| item.1 == variant.1).count();
                if new_count >= self.max_calculate_count || type_count >= self.partition_size {
                    false
                } else {
                    let close = partition.iter().all(|item| {
                        get_variant(variant.0, variant.1)
                            .zip(get_variant(item.0, item.1))
                            .is_some_and(|(current, previous)| {
                                current.start - previous.end < self.max_candidates_distance
                            })
                    });
                    let last_end = partition
                        .last()
                        .and_then(|item| get_variant(item.0, item.1))
                        .map_or(0, |variant| variant.end);
                    close && variant.0 - last_end < EXTEND_BP
                }
            };
            if matches {
                partition.push(variant);
                product = new_count;
                continue;
            }
            if partition
                .last()
                .is_some_and(|previous| previous.0 == variant.0 && previous.1 != variant.1)
            {
                partition.push(variant);
                if let Some(next) = variants.get(index + 1)
                    && !truth_positions_set.contains(&next.0)
                    && next.0 - variant.0 <= EXTEND_BP
                {
                    partition.push(*next);
                }
                partitions.push(std::mem::take(&mut partition));
                product = 1;
            } else {
                partitions.push(std::mem::take(&mut partition));
                partition.push(variant);
                product = combination_count(variant.0, variant.1);
            }
        }
        if !partition.is_empty() {
            partitions.push(partition);
        }
        partitions
            .into_iter()
            .map(|partition| {
                let mut candidates = Vec::new();
                let mut truths = Vec::new();
                for (position, kind) in partition {
                    if let Some(variant) = get_variant(position, kind) {
                        if kind == VariantKind::Candidate {
                            candidates.push(variant);
                        } else {
                            truths.push(variant);
                        }
                    }
                }
                (candidates, truths)
            })
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    fn unify_regions(
        &self,
        candidates: &[Position],
        truths: &[Position],
        split: (i64, i64),
        requested: (Option<i64>, Option<i64>),
        alt_positions: &BTreeMap<i64, Position>,
        truth_positions: &BTreeMap<i64, Position>,
        reads: &BTreeMap<String, ReadInfo>,
        output: &mut String,
        all_positions: &mut BTreeSet<i64>,
        rescued: &mut BTreeMap<i64, String>,
    ) -> Result<()> {
        let record = |position: &Position, genotype: Genotype| {
            format!(
                "{}\t{}\t.\t{}\t{}\t10\tPASS\t.\tGT:GQ:DP:AF\t{}/{}:10:10:0.5000",
                self.contig,
                position.start,
                position.reference,
                position.alternates.join(","),
                genotype.0,
                genotype.1
            )
        };
        for (candidates, truths) in
            self.partitions(candidates, truths, alt_positions, truth_positions)
        {
            let reference = self.reference_window(&candidates, &truths)?;
            let matched =
                self.find_match_pairs(&candidates, &truths, &reference, truth_positions, reads);
            let in_split = |position: i64| {
                position >= split.0
                    && position < split.1
                    && match requested {
                        (Some(start), Some(end)) => position >= start && position < end,
                        _ => true,
                    }
            };
            let Some(matched) = matched else {
                for truth in truths {
                    if !in_split(truth.start) {
                        continue;
                    }
                    let Some(candidate) = alt_positions.get(&truth.start) else {
                        continue;
                    };
                    if !match_alt_base(candidate, &truth)
                        || all_positions.contains(&truth.start)
                        || rescued.contains_key(&truth.start)
                    {
                        continue;
                    }
                    rescued.insert(truth.start, record(&truth, truth.genotype));
                }
                continue;
            };
            let missing_truth = matched
                .truth_genotypes
                .iter()
                .any(|genotype| genotype.0 + genotype.1 == 0);
            for (candidate, genotype) in matched.candidates.iter().zip(&matched.candidate_genotypes)
            {
                if genotype.0 + genotype.1 == 0
                    && missing_truth
                    && !truth_positions.contains_key(&candidate.start)
                    && let Some(phased) = candidate.phased
                    && !rescued.contains_key(&candidate.start)
                {
                    rescued.insert(candidate.start, record(candidate, phased));
                    continue;
                }
                if genotype.0 + genotype.1 == 0 || candidate.alternates.is_empty() {
                    continue;
                }
                let mut selected = BTreeSet::new();
                for allele in [genotype.0, genotype.1] {
                    if allele != 0
                        && let Some(alternate) = candidate.alternates.get(allele as usize - 1)
                    {
                        selected.insert(alternate.clone());
                    }
                }
                let genotype_string = if genotype.0 == 0 || genotype.1 == 0 {
                    (0, 1)
                } else if genotype.0 == genotype.1 {
                    (1, 1)
                } else {
                    (1, 2)
                };
                let (reference, alternates) = remove_common_suffix(
                    &candidate.reference,
                    &selected.into_iter().collect::<Vec<_>>(),
                );
                if alternates.is_empty() {
                    continue;
                }
                let normalized =
                    Position::new(candidate.start, reference, alternates, genotype_string);
                if all_positions.insert(candidate.start) {
                    writeln!(output, "{}", record(&normalized, genotype_string))?;
                    rescued.remove(&candidate.start);
                }
            }
            for (truth, genotype) in matched.truths.iter().zip(&matched.truth_genotypes) {
                if !in_split(truth.start)
                    || *genotype != (0, 0)
                    || truth.genotype.0 + truth.genotype.1 == 0
                {
                    continue;
                }
                let Some(candidate) = alt_positions.get(&truth.start) else {
                    continue;
                };
                if match_alt_base(candidate, truth) && all_positions.insert(truth.start) {
                    rescued.insert(truth.start, record(truth, truth.genotype));
                }
            }
        }
        Ok(())
    }
}

pub fn candidate_details_from_bam(config: &CallConfig, reference: &Reference) -> Result<String> {
    let region = calling_region(reference, config)?;
    if region.start == region.end {
        return Ok(String::new());
    }
    let confident_regions = config
        .extend_bed
        .as_ref()
        .or(config.bed.as_ref())
        .map(|path| IntervalMap::from_bed(path, Some(&config.contig)))
        .transpose()?;
    let known_positions = config
        .known_vcf
        .as_ref()
        .map(|path| candidate_positions(path, Some(&config.contig)))
        .transpose()?
        .map(|positions| positions.into_iter().collect::<BTreeSet<_>>());
    let candidates = create_pileup(
        &config.bam,
        reference,
        &region,
        &PileupOptions {
            min_depth: config.min_coverage,
            min_snp_af: config.min_af.unwrap_or(config.snp_min_af),
            min_indel_af: config.min_af.unwrap_or(config.indel_min_af),
            min_mq: config.min_mq,
            allow_empty_flanking: true,
            call_head_tail: true,
            confident_regions: confident_regions.as_ref(),
            known_positions: known_positions.as_ref(),
            ..Default::default()
        },
    )?;
    let positions = candidates
        .positions
        .iter()
        .map(|position| {
            position
                .rsplit(':')
                .nth(1)
                .context("invalid candidate position")?
                .parse::<i64>()
                .map(|position| position - 1)
                .map_err(anyhow::Error::from)
        })
        .collect::<Result<Vec<_>>>()?;
    create_candidate_details(
        &config.bam,
        reference,
        &config.contig,
        &positions,
        config.platform == Platform::Ilmn,
    )
}

fn candidate_rows(args: &UnifyRepresentationArgs) -> Result<Vec<CandidateRow>> {
    let mut lines = Vec::new();
    if let Some(path) = &args.candidate_details {
        if path.exists() {
            lines.extend(read_text_lines(path)?.collect::<Result<Vec<_>, _>>()?);
        } else {
            let directory = path.parent().unwrap_or_else(|| Path::new("."));
            let prefix = path
                .file_name()
                .context("candidate detail prefix has no filename")?
                .to_string_lossy();
            let mut paths: Vec<PathBuf> = fs::read_dir(directory)?
                .map(|entry| entry.map(|entry| entry.path()))
                .collect::<std::io::Result<_>>()?;
            paths.sort();
            for path in paths {
                if path
                    .file_name()
                    .is_some_and(|name| name.to_string_lossy().starts_with(prefix.as_ref()))
                {
                    lines.extend(read_text_lines(&path)?.collect::<Result<Vec<_>, _>>()?);
                }
            }
        }
    } else {
        let reference = Reference::open(&args.reference)?;
        let bam = args
            .bam
            .as_ref()
            .context("--bam_fn is required when --candidate_details_fn is not provided")?;
        let mut config = CallConfig::new(
            bam.clone(),
            args.reference.clone(),
            PathBuf::new(),
            args.contig.clone(),
            args.platform,
        );
        config.start = args.start;
        config.end = args.end;
        config.bed = args.bed.clone();
        config.extend_bed = args.extend_bed.clone();
        config.known_vcf = args.known_vcf.clone();
        config.chunk_id = args.chunk_id;
        config.chunk_num = args.chunk_num;
        config.min_af = args.min_af;
        config.min_coverage = args.min_coverage;
        config.min_mq = args.min_mq;
        let details = candidate_details_from_bam(&config, &reference)?;
        lines.extend(details.lines().map(str::to_string));
    }
    let mut output = Vec::new();
    for line in lines {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            bail!("invalid candidate detail row: {line}");
        }
        let mut location = fields[0].split_whitespace();
        let contig = location.next().context("candidate contig is missing")?;
        let position = location
            .next()
            .context("candidate position is missing")?
            .parse::<i64>()?;
        if contig != args.contig {
            continue;
        }
        output.push((
            position,
            fields[1].parse()?,
            serde_json::from_str(fields[2])?,
        ));
    }
    output.sort_by_key(|(position, _, _)| *position);
    output.dedup_by_key(|(position, _, _)| *position);
    Ok(output)
}

fn vcf_header(sample_name: &str) -> String {
    format!(
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FILTER=<ID=LowQual,Description=\"Confidence in this variant being real is below calling threshold.\">\n\
##FILTER=<ID=RefCall,Description=\"Genotyping model thinks this site is reference.\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n\
##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
    )
}

fn unified_vcf(args: &UnifyRepresentationArgs) -> Result<String> {
    if args.max_calculate_count == 0 || args.partition_size == 0 {
        bail!("partition size and maximum calculation count must be positive");
    }
    if !(0.0..=1.0).contains(&args.minimum_allele_gap) {
        bail!("minimum_allele_gap must be between zero and one");
    }
    if args.start.is_some() != args.end.is_some() {
        bail!("ctgStart and ctgEnd must be supplied together");
    }
    let requested = (args.start, args.end);
    let candidates = candidate_rows(args)?;
    let (mut start, mut end) = requested;
    if (start.is_none() || end.is_none()) && !candidates.is_empty() {
        if args.candidate_details.is_some() {
            let first = candidates.first().unwrap().0;
            let last = candidates.last().unwrap().0;
            let chunk_count = args
                .chunk_num
                .context("--chunk_num is required when candidate details have no contig range")?;
            let chunk_id = args
                .chunk_id
                .context("--chunk_id is required when candidate details have no contig range")?;
            if chunk_id == 0 || chunk_id > chunk_count {
                bail!("chunk ID must be between 1 and the chunk count");
            }
            let size = (last - first) / chunk_count as i64 + 1;
            let chunk_start = first + size * (chunk_id - 1) as i64;
            start = Some(chunk_start);
            end = Some(chunk_start + size);
        } else {
            let reference = Reference::open(&args.reference)?;
            let length = reference.contig_len(&args.contig) as i64;
            start = Some(args.start.unwrap_or(1));
            end = Some(args.end.unwrap_or(length));
        }
    }
    let first = candidates.first().map_or(start.unwrap_or(1), |row| row.0);
    let last = candidates.last().map_or(end.unwrap_or(first), |row| row.0);
    let reference_start = (start.unwrap_or(first) - REFERENCE_REGION_SIZE).max(1);
    let reference = Reference::open(&args.reference)?;
    let reference_end = (end.unwrap_or(last) + REFERENCE_REGION_SIZE)
        .min(reference.contig_len(&args.contig) as i64);
    let reference_sequence = reference.fetch(&args.contig, reference_start - 1, reference_end)?;
    let mut truth_positions = BTreeMap::new();
    for line in read_text_lines(&args.truth)? {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 || fields[0] != args.contig {
            continue;
        }
        let position = fields[1].parse::<i64>()?;
        if start.zip(end).is_some_and(|(start, end)| {
            position < start - EXTENDED_WINDOW || position > end + EXTENDED_WINDOW
        }) {
            continue;
        }
        let genotype = (fields[4].parse()?, fields[5].parse()?);
        if fields[2] == "None" || fields[3] == "None" || genotype.0 < 0 || genotype.1 < 0 {
            continue;
        }
        truth_positions.insert(
            position,
            Position::new(
                position,
                fields[2].to_string(),
                fields[3].split(',').map(str::to_string).collect(),
                genotype,
            ),
        );
    }
    let mut alt_positions = BTreeMap::new();
    let mut reads = BTreeMap::<String, ReadInfo>::new();
    for (position, depth, candidate_reads) in &candidates {
        if candidate_reads.is_empty()
            || start.zip(end).is_some_and(|(start, end)| {
                *position < start - EXTENDED_WINDOW || *position > end + EXTENDED_WINDOW
            })
        {
            continue;
        }
        let reference_base = reference_sequence
            .get((*position - reference_start) as usize)
            .copied()
            .context("candidate position is outside the reference")?
            as char;
        let (reference, alternates, alt_list) = decode_alt_info(
            candidate_reads,
            &reference_base.to_string(),
            *depth,
            args.minimum_allele_gap,
        );
        let mut candidate = Position::new(*position, reference, alternates, (-1, -1));
        candidate.depth = *depth;
        candidate.alt_list = alt_list;
        for (allele, read_list) in candidate_reads {
            for read in read_list.split(' ').filter(|read| !read.is_empty()) {
                let (name, haplotype) = read
                    .rsplit_once('_')
                    .with_context(|| format!("candidate read {read} has no haplotype suffix"))?;
                let haplotype = haplotype.parse::<usize>()?;
                let read_info = reads.entry(name.to_string()).or_default();
                if read_info.haplotype == 0 && haplotype != 0 {
                    read_info.haplotype = haplotype;
                }
                let haplotype = read_info.haplotype;
                *candidate
                    .variant_haplotypes
                    .entry(allele.clone())
                    .or_default()
                    .entry(haplotype)
                    .or_default() += 1;
                *candidate.haplotype_counts.entry(haplotype).or_default() += 1;
                candidate.reads.insert(name.to_string());
                read_info.position_alleles.insert(*position, allele.clone());
            }
        }
        let truth = truth_positions.get(position);
        let confident = lock_variant(&candidate, truth);
        if confident && let Some(truth) = truth_positions.get_mut(position) {
            truth.confident = true;
        }
        candidate.phased = phased_genotype(
            &candidate,
            truth_positions.get(position),
            reference_base as u8,
            args.minimum_allele_gap,
            args.platform,
        );
        if truth_positions.contains_key(position) && confident {
            if let Some(phased) = candidate.phased {
                let truth = truth_positions.get(position).unwrap();
                let truth_zero_count = [truth.genotype.0, truth.genotype.1]
                    .iter()
                    .filter(|allele| **allele == 0)
                    .count();
                let phased_zero_count = [phased.0, phased.1]
                    .iter()
                    .filter(|allele| **allele == 0)
                    .count();
                if truth_zero_count != phased_zero_count
                    || (truth.genotype.0 + truth.genotype.1 == 3 && phased.0 + phased.1 != 3)
                {
                    candidate.phased = None;
                }
            }
            let truth = truth_positions.get_mut(position).unwrap();
            truth.set_representation(candidate.reference.clone(), candidate.alternates.clone());
            truth.phased = candidate.phased;
        }
        alt_positions.insert(*position, candidate);
    }
    for read in reads.values_mut() {
        for (position, allele) in &read.position_alleles {
            match allele.as_bytes()[0] {
                b'X' => read
                    .sequences
                    .push((*position, *position + 1, allele[1..].to_string())),
                b'I' => read
                    .sequences
                    .push((*position, *position + 1, allele[1..].to_string())),
                b'D' => {
                    let reference_base = reference_sequence
                        .get((*position - reference_start) as usize)
                        .copied()
                        .context("deletion position is outside the reference")?;
                    read.sequences.push((
                        *position,
                        *position + allele.len() as i64,
                        (reference_base as char).to_string(),
                    ));
                }
                _ => {}
            }
        }
    }
    if alt_positions.is_empty() || truth_positions.is_empty() {
        return Ok(vcf_header(&args.sample_name));
    }
    let confident_regions = args
        .bed
        .as_ref()
        .map(|path| IntervalMap::from_bed(path, Some(&args.contig)))
        .transpose()?;
    let region_start = start.unwrap_or_else(|| alt_positions.keys().min().copied().unwrap());
    let region_end = end.unwrap_or_else(|| alt_positions.keys().max().copied().unwrap());
    let unification = Unification {
        contig: &args.contig,
        reference: &reference_sequence,
        reference_start,
        partition_size: args.partition_size,
        max_candidates_distance: if args.platform == Platform::Hifi {
            200
        } else {
            args.max_candidates_distance
        },
        max_calculate_count: args.max_calculate_count,
    };
    let mut partition_size = args.partition_size;
    if args.platform == Platform::Hifi {
        partition_size = 20;
    }
    let unification = Unification {
        partition_size,
        ..unification
    };
    let _ = (args.subsample_ratio, args.test_pos, &args.pypy);
    let mut output = String::new();
    let mut all_positions = BTreeSet::new();
    let mut rescued = BTreeMap::new();
    for index in 0..=(region_end - region_start) / REGION_SIZE {
        let split_start = region_start + index * REGION_SIZE;
        let split_end = split_start + REGION_SIZE;
        let candidates: Vec<Position> = alt_positions
            .values()
            .filter(|candidate| {
                candidate.start >= split_start - EXTEND_BP
                    && candidate.start < split_end + EXTEND_BP
                    && !candidate.alternates.is_empty()
                    && confident_regions.as_ref().is_none_or(|regions| {
                        regions.overlaps(&args.contig, candidate.start - 2, candidate.end + 2)
                    })
            })
            .cloned()
            .collect();
        let truths: Vec<Position> = truth_positions
            .values()
            .filter(|truth| {
                truth.start >= split_start - EXTEND_BP
                    && truth.start < split_end + EXTEND_BP
                    && confident_regions.as_ref().is_none_or(|regions| {
                        regions.overlaps(&args.contig, truth.start - 2, truth.end + 2)
                    })
            })
            .cloned()
            .collect();
        if candidates.is_empty() && truths.is_empty() {
            continue;
        }
        unification.unify_regions(
            &candidates,
            &truths,
            (split_start, split_end),
            requested,
            &alt_positions,
            &truth_positions,
            &reads,
            &mut output,
            &mut all_positions,
            &mut rescued,
        )?;
    }
    for record in rescued.values() {
        writeln!(output, "{record}")?;
    }
    Ok(format!("{}{}", vcf_header(&args.sample_name), output))
}

pub fn unify_representation(args: UnifyRepresentationArgs) -> Result<()> {
    let output = unified_vcf(&args)?;
    if let Some(path) = &args.output {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, output)?;
    } else {
        print!("{output}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn suffix_trimming_and_genotypes_match_unification_rules() {
        assert_eq!(
            remove_common_suffix("ACGT", &["AGT".to_string(), "ACGGT".to_string()]),
            ("AC".to_string(), vec!["A".to_string(), "ACG".to_string()])
        );
        let options = vec![
            BTreeSet::from([(0, 0), (0, 1), (1, 0)]),
            BTreeSet::from([(0, 0), (0, 1), (1, 0)]),
        ];
        let combinations = unique_genotypes(&options);
        assert_eq!(combinations.len(), 5);
        assert!(combinations.iter().all(|genotypes| genotypes.len() == 2));
    }

    #[test]
    fn haplotypes_distinguish_shifted_insertions() {
        let reference = ReferenceWindow {
            start: 4,
            end: 16,
            sequence: b"AACCCCCCCCCC".to_vec(),
        };
        let first = Position::new(5, "A".to_string(), vec!["AC".to_string()], (0, 1));
        let second = Position::new(6, "C".to_string(), vec!["CC".to_string()], (0, 1));
        let first_sequences = sequence_combinations(
            std::slice::from_ref(&first),
            vec![vec![(0, 1)], vec![(1, 0)]],
            &reference,
        );
        let second_sequences = sequence_combinations(
            std::slice::from_ref(&second),
            vec![vec![(0, 1)], vec![(1, 0)]],
            &reference,
        );
        assert_eq!(first_sequences.len(), 1);
        assert!(first_sequences.keys().all(|sequences| sequences.len() == 2));
        let homozygous_sequences =
            sequence_combinations(std::slice::from_ref(&first), vec![vec![(1, 1)]], &reference);
        assert!(
            homozygous_sequences
                .keys()
                .all(|sequences| sequences.len() == 1)
        );
        assert_eq!(first_sequences, second_sequences);
    }

    #[test]
    fn bam_candidate_discovery_includes_long_indels() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data");
        let reference = Reference::open(&data.join("toy.fa")).unwrap();
        let mut config = CallConfig::new(
            data.join("toy.bam"),
            data.join("toy.fa"),
            PathBuf::new(),
            "chrToy".to_string(),
            Platform::Ont,
        );
        config.start = Some(1);
        config.end = Some(4096);
        let details = candidate_details_from_bam(&config, &reference).unwrap();
        let positions: Vec<i64> = details
            .lines()
            .map(|line| line.split_whitespace().nth(1).unwrap().parse().unwrap())
            .collect();
        assert_eq!(
            positions,
            [
                5, 500, 700, 900, 1100, 1500, 1900, 2300, 2500, 2700, 3100, 4090
            ]
        );
    }

    fn test_arguments(
        reference: PathBuf,
        truth: PathBuf,
        candidates: PathBuf,
        bed: Option<PathBuf>,
        end: i64,
    ) -> UnifyRepresentationArgs {
        UnifyRepresentationArgs {
            platform: Platform::Ont,
            truth,
            reference,
            candidate_details: Some(candidates),
            output: None,
            sample_name: "SAMPLE".to_string(),
            contig: "chrToy".to_string(),
            start: Some(1),
            end: Some(end),
            max_candidates_distance: 100,
            max_calculate_count: 10_000,
            partition_size: 15,
            minimum_allele_gap: 0.15,
            bed,
            known_vcf: None,
            extend_bed: None,
            subsample_ratio: 1000,
            test_pos: 0,
            chunk_num: None,
            chunk_id: None,
            min_af: None,
            bam: None,
            pypy: String::new(),
            min_coverage: 2,
            min_mq: 5,
        }
    }

    #[test]
    fn direct_bam_unification_matches_toy_truth() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data");
        let reference = Reference::open(&data.join("toy.fa")).unwrap();
        let mut config = CallConfig::new(
            data.join("toy.bam"),
            data.join("toy.fa"),
            PathBuf::new(),
            "chrToy".to_string(),
            Platform::Ont,
        );
        config.start = Some(1);
        config.end = Some(4096);
        config.bed = Some(data.join("small.bed"));
        let details = candidate_details_from_bam(&config, &reference).unwrap();
        let directory = tempfile::tempdir().unwrap();
        let candidates = directory.path().join("candidates.txt");
        fs::write(&candidates, details).unwrap();
        let output = unified_vcf(&test_arguments(
            data.join("toy.fa"),
            data.join("truth.var.txt"),
            candidates,
            Some(data.join("small.bed")),
            4096,
        ))
        .unwrap();
        let records: Vec<(i64, String, String, String)> = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                (
                    fields[1].parse().unwrap(),
                    fields[3].to_string(),
                    fields[4].to_string(),
                    fields[9].split(':').next().unwrap().to_string(),
                )
            })
            .collect();
        assert_eq!(
            records,
            [
                (500, "A".to_string(), "C".to_string(), "0/1".to_string()),
                (700, "T".to_string(), "TGGG".to_string(), "0/1".to_string()),
                (900, "GCA".to_string(), "G".to_string(), "0/1".to_string()),
                (2500, "T".to_string(), "A".to_string(), "0/1".to_string()),
                (1100, "A".to_string(), "C,G".to_string(), "1/2".to_string()),
            ]
        );
    }

    #[test]
    fn homopolymer_ties_match_the_legacy_nonvariant_tie_break() {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chrToy\n{}\n", "A".repeat(200))).unwrap();
        fs::write(
            reference.with_extension("fa.fai"),
            "chrToy\t200\t8\t200\t201\n",
        )
        .unwrap();
        let truth = directory.path().join("truth.var");
        fs::write(&truth, "chrToy 100 A AA 1 1\n").unwrap();
        let candidates = directory.path().join("candidates.txt");
        fs::write(
            &candidates,
            concat!(
                "chrToy 100\t2\t{\"IAA\": \"low1_1 low2_2\"}\n",
                "chrToy 101\t4\t{\"IAA\": \"high1_1 high2_1 high3_2 high4_2\"}\n"
            ),
        )
        .unwrap();
        let bed = directory.path().join("confident.bed");
        fs::write(&bed, "chrToy\t0\t200\n").unwrap();
        let output = unified_vcf(&test_arguments(
            reference,
            truth,
            candidates,
            Some(bed),
            200,
        ))
        .unwrap();
        let records: Vec<&str> = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect();
        assert_eq!(records.len(), 2, "{output}");
        for (index, position) in [100, 101].into_iter().enumerate() {
            assert!(
                records[index].starts_with(&format!("chrToy\t{position}\t.\tA\tAA\t")),
                "{output}"
            );
            assert!(records[index].ends_with("0/1:10:10:0.5000"), "{output}");
        }
    }
}
