use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

use crate::{
    io::{reference::Reference, variant::VcfRecord},
    pipeline::merge::read_vcf,
};
use anyhow::{Context, Result, bail};

#[derive(Debug, Clone)]
pub struct SelectOptions {
    pub var_pct_full: f64,
    pub ref_pct_full: Option<f64>,
    pub var_pct_phasing: f64,
    pub split_bed_size: usize,
    pub phasing_window_size: i64,
    pub phasing_info_in_bam: bool,
    pub call_low_seq_entropy: bool,
    pub seq_entropy_pro: f64,
    pub qual_path: Option<PathBuf>,
    pub all_alt_path: Option<PathBuf>,
}

impl Default for SelectOptions {
    fn default() -> Self {
        Self {
            var_pct_full: 0.3,
            ref_pct_full: None,
            var_pct_phasing: 0.7,
            split_bed_size: 10000,
            phasing_window_size: 100000,
            phasing_info_in_bam: true,
            call_low_seq_entropy: false,
            seq_entropy_pro: 0.05,
            qual_path: None,
            all_alt_path: None,
        }
    }
}

fn cutoff(mut values: Vec<f64>, proportion: f64) -> f64 {
    values.sort_by(|left, right| left.total_cmp(right));
    let count = (proportion * values.len() as f64) as usize;
    if count == 0 { 0.0 } else { values[count - 1] }
}

pub fn phase_quality(
    records: &[VcfRecord],
    var_pct_phasing: Option<f64>,
    var_pct_full: f64,
) -> f64 {
    let proportion = var_pct_phasing.map_or(var_pct_full, |value| 1.0 - value);
    let qualities: Vec<f64> = records
        .iter()
        .filter(|record| {
            record.reference.len() == 1
                && record.alternate.len() == 1
                && record.samples[0].split(':').next() == Some("0/1")
        })
        .filter_map(|record| record.quality.parse().ok())
        .collect();
    cutoff(qualities, proportion)
}

pub fn full_quality(
    records: &[VcfRecord],
    var_pct_full: f64,
    ref_pct_full: Option<f64>,
    known_vcf: bool,
) -> (f64, f64) {
    let variants: Vec<f64> = records
        .iter()
        .filter(|record| record.samples[0].split(':').next() != Some("0/0"))
        .filter_map(|record| record.quality.parse().ok())
        .collect();
    let references: Vec<f64> = records
        .iter()
        .filter(|record| record.samples[0].split(':').next() == Some("0/0"))
        .filter_map(|record| record.quality.parse().ok())
        .collect();
    let reference_proportion = if known_vcf {
        ref_pct_full.unwrap_or(0.0)
    } else {
        ref_pct_full
            .filter(|value| *value != 0.0)
            .unwrap_or(var_pct_full)
    };
    (
        cutoff(variants, var_pct_full),
        cutoff(references, reference_proportion),
    )
}

pub fn select_het_snp(
    records: &[VcfRecord],
    contig: &str,
    var_pct_full: f64,
    phase_qual: Option<f64>,
) -> Vec<VcfRecord> {
    let mut candidates = BTreeMap::new();
    let mut qualities = BTreeMap::new();
    for record in records {
        if record.contig != contig
            || record.reference.len() != 1
            || record.alternate.len() != 1
            || !matches!(record.samples[0].split(':').next(), Some("0/1" | "1/0"))
        {
            continue;
        }
        let mut record = record.clone();
        if let Some(index) = record
            .format
            .as_ref()
            .and_then(|format| format.split(':').position(|field| field == "PL"))
        {
            let format_string = record.format.unwrap();
            let mut format: Vec<_> = format_string.split(':').collect();
            format.remove(index);
            record.format = Some(format.join(":"));
            let mut sample: Vec<_> = record.samples[0].split(':').collect();
            sample.remove(index);
            record.samples[0] = sample.join(":");
        }
        let position = record.position;
        let quality = record.quality.parse().unwrap_or(0.0);
        candidates.insert(position, record);
        qualities.insert(position, quality);
    }
    let remove: BTreeSet<i64> = if let Some(threshold) = phase_qual {
        qualities
            .iter()
            .filter(|(_, quality)| **quality < threshold)
            .map(|(position, _)| *position)
            .collect()
    } else {
        let mut values: Vec<(i64, f64)> = qualities.into_iter().collect();
        values.sort_by(|left, right| left.1.partial_cmp(&right.1).unwrap());
        let count = (var_pct_full * values.len() as f64) as usize;
        values
            .into_iter()
            .take(count)
            .map(|(position, _)| position)
            .collect()
    };
    candidates
        .into_iter()
        .filter(|(position, _)| !remove.contains(position))
        .map(|(_, record)| record)
        .collect()
}

fn base_number(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' | b'R' | b'W' | b'M' | b'D' | b'H' | b'V' | b'N' => 0,
        b'C' | b'Y' | b'S' | b'B' => 1,
        b'G' | b'K' => 2,
        b'T' | b'U' => 3,
        _ => 0,
    }
}

fn sequence_entropy(sequence: &[u8], entropy_window: usize, kmer: usize) -> f64 {
    let mut repeat_counts = vec![0_i64; entropy_window + 2];
    repeat_counts[0] = entropy_window as i64;
    let mut entropy = vec![0.0; entropy_window + 2];
    for (index, value) in entropy.iter_mut().enumerate().skip(1) {
        let e = index as f64 / entropy_window as f64;
        *value = e * e.ln();
    }
    let entropy_mul = -1.0 / (entropy_window as f64).ln();
    let mask = (1_usize << (2 * kmer)) - 1;
    let mut hash_counts = vec![0_i64; 1 << (2 * kmer)];
    let mut suffix = 0_usize;
    let mut prefix = 0_usize;
    let mut entropy_sum = 0.0;
    let mut index = 0_isize;
    let mut previous = -(entropy_window as isize);
    while previous < sequence.len() as isize {
        if index < sequence.len() as isize {
            suffix = ((suffix << 2) | base_number(sequence[index as usize])) & mask;
            let old = hash_counts[suffix] as usize;
            repeat_counts[old] -= 1;
            entropy_sum -= entropy[old];
            hash_counts[suffix] += 1;
            let new = hash_counts[suffix] as usize;
            repeat_counts[new] += 1;
            entropy_sum += entropy[new];
        }
        if previous >= 0 && index < sequence.len() as isize {
            prefix = ((prefix << 2) | base_number(sequence[previous as usize])) & mask;
            let old = hash_counts[prefix] as usize;
            repeat_counts[old] -= 1;
            entropy_sum -= entropy[old];
            hash_counts[prefix] -= 1;
            let new = hash_counts[prefix] as usize;
            repeat_counts[new] += 1;
            entropy_sum += entropy[new];
        }
        index += 1;
        previous += 1;
    }
    entropy_sum * entropy_mul
}

fn entropy_positions(
    reference: &Reference,
    contig: &str,
    positions: &[i64],
    proportion: f64,
) -> Result<Vec<i64>> {
    if positions.is_empty() {
        return Ok(Vec::new());
    }
    let start = (*positions.iter().min().unwrap() - 33).max(1);
    let end = *positions.iter().max().unwrap() + 34;
    let sequence = reference.fetch(contig, start - 1, end)?;
    let mut entropy = Vec::new();
    for position in positions {
        let sequence =
            &sequence[(position - 16 - start) as usize..(position + 17 - start) as usize];
        entropy.push((*position, sequence_entropy(sequence, 33, 5)));
    }
    entropy.sort_by(|left, right| left.1.partial_cmp(&right.1).unwrap());
    let count = (proportion * entropy.len() as f64) as usize;
    Ok(entropy
        .into_iter()
        .take(count)
        .map(|(position, _)| position)
        .collect())
}

fn write_regions(output_dir: &Path, contig: &str, files: &[PathBuf]) -> Result<()> {
    let regions = output_dir.join(format!("full_aln_regions_{contig}"));
    let mut content = String::new();
    for path in files {
        content.push_str(&fs::read_to_string(path)?);
    }
    fs::write(regions, content)?;
    fs::write(
        output_dir.join(format!("FULL_ALN_FILE_{contig}")),
        files
            .iter()
            .map(|path| path.display().to_string())
            .collect::<Vec<_>>()
            .join("\n")
            + if files.is_empty() { "" } else { "\n" },
    )?;
    Ok(())
}

pub fn select_candidates(
    pileup_vcf: &Path,
    phased_vcf: Option<&Path>,
    output_dir: &Path,
    contig: &str,
    platform: &str,
    options: &SelectOptions,
    reference: Option<&Reference>,
) -> Result<Vec<PathBuf>> {
    if options.split_bed_size == 0 {
        bail!("split BED size must be positive");
    }
    fs::create_dir_all(output_dir)?;
    let mut references = Vec::new();
    let mut variants = Vec::new();
    let mut positions = BTreeSet::new();
    let has_pileup = pileup_vcf.exists();
    if has_pileup {
        for record in &read_vcf(pileup_vcf)?.records {
            if record.contig != contig {
                continue;
            }
            let quality = record.quality.parse::<f64>()?;
            if record.alternate == "." || record.reference == record.alternate {
                references.push((record.position, quality));
            } else {
                variants.push((record.position, quality));
                positions.insert(record.position);
            }
        }
    } else if let Some(path) = &options.all_alt_path {
        for line in crate::io::compression::read_text_lines(path)? {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() < 2 || fields[0] != contig {
                continue;
            }
            positions.insert(fields[1].parse()?);
        }
    }
    let reference_proportion = options.ref_pct_full.unwrap_or(options.var_pct_full);
    let (var_cutoff, ref_cutoff) = if let Some(path) = &options.qual_path {
        let text = fs::read_to_string(path)?;
        let values: Vec<f64> = text
            .split_whitespace()
            .map(|value| value.parse())
            .collect::<Result<_, _>>()?;
        if values.len() != 2 {
            bail!("quality cutoff file must contain variant and reference cutoffs");
        }
        (values[0], values[1])
    } else {
        (
            cutoff(
                variants.iter().map(|(_, quality)| *quality).collect(),
                options.var_pct_full,
            ),
            cutoff(
                references.iter().map(|(_, quality)| *quality).collect(),
                reference_proportion,
            ),
        )
    };
    let low_variants: Vec<i64> = variants
        .iter()
        .filter(|(_, quality)| *quality < var_cutoff)
        .map(|(position, _)| *position)
        .collect();
    let low_references: Vec<i64> = references
        .iter()
        .filter(|(_, quality)| *quality < ref_cutoff)
        .map(|(position, _)| *position)
        .collect();
    let mut low_entropy = Vec::new();
    if options.call_low_seq_entropy && has_pileup {
        let proportion = options.var_pct_full + options.seq_entropy_pro;
        let mut candidates: Vec<(i64, f64)> = references
            .iter()
            .copied()
            .sorted_by_quality()
            .take((proportion * references.len() as f64) as usize)
            .collect();
        candidates.extend(
            variants
                .iter()
                .copied()
                .sorted_by_quality()
                .take((proportion * variants.len() as f64) as usize),
        );
        candidates.sort_by_key(|(position, _)| *position);
        let candidate_positions: Vec<i64> =
            candidates.iter().map(|(position, _)| *position).collect();
        if !candidate_positions.is_empty() {
            let reference = reference.context("sequence entropy selection requires a reference")?;
            low_entropy = entropy_positions(
                reference,
                contig,
                &candidate_positions,
                options.seq_entropy_pro,
            )?;
        }
    }
    if options.phasing_info_in_bam && has_pileup {
        positions = low_references.into_iter().collect();
        positions.extend(low_variants);
        positions.extend(low_entropy);
    } else {
        positions.extend(low_references);
        positions.extend(low_entropy);
    }
    let phased = if let Some(path) = phased_vcf.filter(|path| path.exists()) {
        let mut variants = BTreeMap::new();
        for record in read_vcf(path)?.records {
            if record.contig != contig {
                continue;
            }
            let fields: Vec<_> = record.samples[0].split(':').collect();
            let genotype = fields[0];
            if !genotype.contains('|') {
                continue;
            }
            variants.insert(
                record.position,
                format!(
                    "{}-{}-{}-{}",
                    record.reference,
                    record.alternate,
                    if genotype == "0|1" { "1" } else { "2" },
                    fields.last().unwrap()
                ),
            );
        }
        variants
    } else {
        BTreeMap::new()
    };
    let positions: Vec<i64> = positions.into_iter().collect();
    let mut files = Vec::new();
    for (chunk_index, chunk) in positions.chunks(options.split_bed_size).enumerate() {
        let path = if options.phasing_info_in_bam && has_pileup {
            output_dir.join(format!(
                "{contig}.{chunk_index}_{}",
                positions.len().div_ceil(options.split_bed_size)
            ))
        } else {
            output_dir.join(format!("{contig}.{}_{}", chunk[0], chunk[chunk.len() - 1]))
        };
        let start = chunk[0] - options.phasing_window_size;
        let end = chunk[chunk.len() - 1] + options.phasing_window_size;
        let mut rows = Vec::new();
        if !options.phasing_info_in_bam || !has_pileup {
            for (position, metadata) in &phased {
                if *position >= start && *position < end && positions.contains(position) {
                    rows.push(format!(
                        "{contig}\t{}\t{}\t{metadata}",
                        position - 12,
                        position + 10
                    ));
                }
            }
        }
        for position in chunk {
            if platform == "ilmn" && options.phasing_info_in_bam && has_pileup {
                let start = (position / 1000 * 1000 - 33).max(1);
                let end = position / 1000 * 1000 + 1033;
                rows.push(format!("{contig}\t{}\t{}", start - 1, end - 1));
            } else if options.phasing_info_in_bam && has_pileup {
                rows.push(format!("{contig}\t{}\t{}", position - 17, position + 17));
            } else {
                rows.push(format!("{contig}\t{}\t{}", position - 17, position + 16));
            }
        }
        rows.sort_by_key(|row| row.split('\t').nth(1).unwrap().parse::<i64>().unwrap());
        fs::write(&path, rows.join("\n") + "\n")?;
        files.push(path);
    }
    write_regions(output_dir, contig, &files)?;
    Ok(files)
}

#[derive(Debug, Clone)]
pub struct TrainingSelectOptions {
    pub var_pct_full: f64,
    pub ref_pct_full: Option<f64>,
    pub split_bed_size: usize,
    pub phasing_info_in_bam: bool,
    pub realign_window_size: Option<i64>,
    pub chunk_num: Option<usize>,
    pub sample_name: String,
    pub depth: usize,
    pub chr_prefix: String,
    pub phased_vcf: Option<PathBuf>,
    pub all_alt_path: Option<PathBuf>,
}

fn write_training_chunk(
    path: &Path,
    positions: &[i64],
    contig: &str,
    options: &TrainingSelectOptions,
    phased: &BTreeMap<i64, String>,
    selected: &BTreeSet<i64>,
) -> Result<()> {
    let window = options
        .realign_window_size
        .unwrap_or(crate::config::FLANKING_BASES as i64);
    let mut rows = Vec::new();
    for position in positions {
        let start = if options.phasing_info_in_bam {
            (position - window - 1).max(0)
        } else {
            position - window - 1
        };
        let end = position + window + i64::from(options.phasing_info_in_bam);
        rows.push((start, format!("{contig}\t{start}\t{end}")));
    }
    if !options.phasing_info_in_bam && !positions.is_empty() {
        let start = positions[0] - 100000;
        let end = positions[positions.len() - 1] + 100000;
        for (position, metadata) in phased.range(start..=end) {
            if selected.contains(position) {
                rows.push((
                    position - 12,
                    format!(
                        "{contig}\t{}\t{}\t{}",
                        position - 12,
                        position + 10,
                        metadata
                    ),
                ));
            }
        }
    }
    rows.sort_by_key(|(start, _)| *start);
    let content = rows
        .into_iter()
        .map(|(_, row)| row)
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    fs::write(path, content)?;
    Ok(())
}

pub fn select_training_candidates(
    alternative_vcf: Option<&Path>,
    output_dir: &Path,
    contig: &str,
    options: &TrainingSelectOptions,
) -> Result<Vec<PathBuf>> {
    if options.split_bed_size == 0 || options.chunk_num == Some(0) {
        bail!("training chunk counts and BED size must be positive");
    }
    fs::create_dir_all(output_dir)?;
    let mut phased = BTreeMap::new();
    if let Some(path) = options.phased_vcf.as_ref().filter(|path| path.exists()) {
        for record in read_vcf(path)?.records {
            if record.contig != contig {
                continue;
            }
            let sample = record.samples.first().context("phased VCF has no sample")?;
            let fields: Vec<_> = sample.split(':').collect();
            let genotype = fields[0];
            if !genotype.contains('|') {
                continue;
            }
            phased.insert(
                record.position,
                format!(
                    "{}-{}-{}-{}",
                    record.reference,
                    record.alternate,
                    if genotype == "0|1" { "1" } else { "2" },
                    fields.last().unwrap()
                ),
            );
        }
    }

    let mut selected = BTreeSet::new();
    if let Some(path) = alternative_vcf {
        let mut references = Vec::new();
        let mut variants = Vec::new();
        for record in read_vcf(path)?.records {
            if record.contig != contig {
                continue;
            }
            let quality = record.quality.parse::<f64>()?;
            if record.alternate == "." || record.reference == record.alternate {
                references.push((record.position, quality));
            } else {
                variants.push((record.position, quality));
            }
        }
        references.sort_by(|left, right| left.1.total_cmp(&right.1));
        variants.sort_by(|left, right| left.1.total_cmp(&right.1));
        let reference_count = (options.ref_pct_full.unwrap_or(options.var_pct_full)
            * references.len() as f64) as usize;
        let variant_count = (options.var_pct_full * variants.len() as f64) as usize;
        let low_references = &references[..reference_count];
        let low_variants = &variants[..variant_count];
        println!(
            "[INFO] {} {} total low qual ref calling to process: {}",
            options.sample_name,
            contig,
            low_references.len()
        );
        println!(
            "[INFO] {} {} total low qual variant calling to process: {}",
            options.sample_name,
            contig,
            low_variants.len()
        );
        if options.phasing_info_in_bam {
            selected.extend(low_references.iter().map(|(position, _)| *position));
            selected.extend(low_variants.iter().map(|(position, _)| *position));
        } else {
            selected.extend(variants.iter().map(|(position, _)| *position));
            selected.extend(low_references.iter().map(|(position, _)| *position));
        }
    } else if let Some(path) = &options.all_alt_path {
        for line in crate::io::compression::read_text_lines(path)? {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() < 2 || fields[0] != contig {
                continue;
            }
            selected.insert(fields[1].parse()?);
        }
    }

    let positions: Vec<i64> = selected.iter().copied().collect();
    println!(
        "[INFO] Total candidates need to be processed in {}: {}",
        contig,
        positions.len()
    );
    let chunks: Vec<&[i64]> = if let Some(chunk_num) = options.chunk_num {
        let chunk_size = positions.len().div_ceil(chunk_num);
        (0..chunk_num)
            .map(|index| {
                let start = index * chunk_size;
                if start < positions.len() {
                    &positions[start..(start + chunk_size).min(positions.len())]
                } else {
                    &[]
                }
            })
            .collect()
    } else {
        positions.chunks(options.split_bed_size).collect()
    };
    let short_contig = contig.strip_prefix(&options.chr_prefix).unwrap_or(contig);
    let mut files = Vec::new();
    for (index, chunk) in chunks.into_iter().enumerate() {
        let path = if options.chunk_num.is_some() {
            output_dir.join(format!(
                "{}_{}_{}_{}",
                options.sample_name,
                options.depth,
                short_contig,
                index + 1
            ))
        } else {
            output_dir.join(format!(
                "{}.{}_{}",
                short_contig,
                chunk[0],
                chunk[chunk.len() - 1]
            ))
        };
        write_training_chunk(&path, chunk, contig, options, &phased, &selected)?;
        files.push(path);
    }
    Ok(files)
}

trait SortByQuality: Iterator<Item = (i64, f64)> + Sized {
    fn sorted_by_quality(self) -> std::vec::IntoIter<(i64, f64)> {
        let mut values: Vec<_> = self.collect();
        values.sort_by(|left, right| left.1.total_cmp(&right.1));
        values.into_iter()
    }
}

impl<T: Iterator<Item = (i64, f64)>> SortByQuality for T {}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> TrainingSelectOptions {
        TrainingSelectOptions {
            var_pct_full: 1.0,
            ref_pct_full: Some(0.5),
            split_bed_size: 1000,
            phasing_info_in_bam: true,
            realign_window_size: None,
            chunk_num: Some(2),
            sample_name: "HG002".to_string(),
            depth: 1000,
            chr_prefix: "chr".to_string(),
            phased_vcf: None,
            all_alt_path: None,
        }
    }

    #[test]
    fn training_selection_matches_chunk_names_and_bed_boundaries() {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("input.vcf");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr21\t100\t.\tA\tC\t20\t.\t.\tGT\t0/1\n\
chr21\t200\t.\tA\tA\t30\t.\t.\tGT\t0/0\n\
chr21\t300\t.\tG\tT\t10\t.\t.\tGT\t1/1\n\
chr21\t400\t.\tC\tC\t10\t.\t.\tGT\t0/0\n\
chr21\t500\t.\tG\tA\t30\t.\t.\tGT\t0/1\n",
        )
        .unwrap();
        let files = select_training_candidates(Some(&input), directory.path(), "chr21", &options())
            .unwrap();
        assert_eq!(files.len(), 2);
        assert_eq!(
            files,
            [
                directory.path().join("HG002_1000_21_1"),
                directory.path().join("HG002_1000_21_2")
            ]
        );
        assert_eq!(
            fs::read_to_string(&files[0]).unwrap(),
            "chr21\t83\t117\nchr21\t283\t317\n"
        );
        assert_eq!(
            fs::read_to_string(&files[1]).unwrap(),
            "chr21\t383\t417\nchr21\t483\t517\n"
        );
    }

    #[test]
    fn training_selection_preserves_phased_metadata() {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("input.vcf");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr21\t100\t.\tA\tC\t20\t.\t.\tGT\t0/1\n",
        )
        .unwrap();
        let phased = directory.path().join("phased.vcf");
        fs::write(
            &phased,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr21\t100\t.\tA\tC\t20\t.\t.\tGT:PS\t0|1:7\n",
        )
        .unwrap();
        let mut options = options();
        options.phasing_info_in_bam = false;
        options.chunk_num = None;
        options.phased_vcf = Some(phased);
        let files =
            select_training_candidates(Some(&input), directory.path(), "chr21", &options).unwrap();
        assert_eq!(files, [directory.path().join("21.100_100")]);
        assert_eq!(
            fs::read_to_string(&files[0]).unwrap(),
            "chr21\t83\t116\nchr21\t88\t110\tA-C-1-7\n"
        );
    }
}
