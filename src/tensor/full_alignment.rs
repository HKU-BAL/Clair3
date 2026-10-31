use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

use super::kh_counter::KhCounter;
use anyhow::Result;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, Record, record::Aux};

use crate::{
    config::{FLANKING_BASES, FULL_ALIGNMENT_CHANNELS, POSITION_COUNT},
    io::{
        compression::read_text_lines, interval::Region, reference::Reference, variant::VcfReader,
    },
    tensor::dwell::signal_lengths,
};

const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const FILTER_FLAGS: u16 = 0x4 | 0x8 | 0x100 | 0x800;
const MIN_HAPLOTAG_MQ: u8 = 20;

struct SamplingRandom(*mut libc::c_void);

impl SamplingRandom {
    fn new(seed: u32) -> Self {
        unsafe extern "C" {
            fn clair3_sampling_new(seed: u32) -> *mut libc::c_void;
        }
        let state = unsafe { clair3_sampling_new(seed) };
        assert!(!state.is_null(), "failed to allocate random state");
        Self(state)
    }

    fn next(&mut self) -> usize {
        unsafe extern "C" {
            fn clair3_sampling_next(state: *mut libc::c_void) -> i32;
        }
        unsafe { clair3_sampling_next(self.0) as usize }
    }
}

impl Drop for SamplingRandom {
    fn drop(&mut self) {
        unsafe extern "C" {
            fn clair3_sampling_free(state: *mut libc::c_void);
        }
        unsafe { clair3_sampling_free(self.0) }
    }
}

#[derive(Debug, Clone)]
pub struct FullAlignmentOptions {
    pub matrix_depth: usize,
    pub max_indel_length: usize,
    pub min_mq: u8,
    pub min_bq: u8,
    pub need_haplotagging: bool,
    pub phasing_info_in_bam: bool,
    pub force_unphased: bool,
    pub enable_dwell_time: bool,
    pub illumina: bool,
}

impl Default for FullAlignmentOptions {
    fn default() -> Self {
        Self {
            matrix_depth: 89,
            max_indel_length: 50,
            min_mq: 5,
            min_bq: 0,
            need_haplotagging: false,
            phasing_info_in_bam: false,
            force_unphased: false,
            enable_dwell_time: false,
            illumina: false,
        }
    }
}

#[derive(Debug, Default)]
pub struct FullAlignmentBatch {
    pub shape: Vec<usize>,
    pub values: Vec<i8>,
    pub positions: Vec<String>,
    pub alt_info: Vec<String>,
}

impl FullAlignmentBatch {
    pub fn append(&mut self, other: Self) -> Result<()> {
        if self.positions.is_empty() {
            *self = other;
            return Ok(());
        }
        if self.shape[1..] != other.shape[1..] {
            anyhow::bail!("cannot append full-alignment tensors with different shapes");
        }
        self.shape[0] += other.shape[0];
        self.values.extend(other.values);
        self.positions.extend(other.positions);
        self.alt_info.extend(other.alt_info);
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct PhasedVariant {
    pub position: i64,
    pub ref_base: u8,
    pub alt_base: u8,
    pub genotype: i32,
    pub phase_set: i32,
}

pub fn phased_variants_from_vcf(path: &Path, contig: &str) -> Result<Vec<PhasedVariant>> {
    let mut reader = VcfReader::open(path)?;
    let mut variants = vec![];
    while let Some(record) = reader.next_record()? {
        if record.contig != contig || record.samples.is_empty() {
            continue;
        }
        let sample = &record.samples[0];
        let fields: Vec<&str> = sample.split(':').collect();
        let genotype = fields[0];
        if !genotype.contains('|') {
            continue;
        }
        let phase_set = fields.last().unwrap().parse()?;
        variants.push(PhasedVariant {
            position: record.position - 1,
            ref_base: record.reference.as_bytes()[0],
            alt_base: record.alternate.as_bytes()[0],
            genotype: if genotype == "0|1" { 1 } else { 2 },
            phase_set,
        });
    }
    Ok(variants)
}

fn candidates_from_regions_with_offset(path: &Path, contig: &str, offset: i64) -> Result<Vec<i64>> {
    let mut candidates = BTreeSet::new();
    let mut contig_start = i64::MAX;
    let mut contig_end = 0_i64;
    for line in read_text_lines(path)? {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.trim_end().split('\t').collect();
        if columns.len() < 3 || columns[0] != contig {
            continue;
        }
        let position: i64 = columns[1].parse::<i64>()? + 1;
        let end: i64 = columns[2].parse::<i64>()? + 1;
        contig_start = contig_start.min(position);
        contig_end = contig_end.max(end);
        if columns.len() == 3 {
            let center = if position == 1 {
                end - FLANKING_BASES as i64 - 2
            } else {
                position + (end - position) / 2 - 1
            };
            candidates.insert(center);
        }
    }
    if candidates.is_empty() {
        return Ok(vec![]);
    }
    Ok(candidates
        .into_iter()
        .filter(|candidate| *candidate >= contig_start && *candidate <= contig_end)
        .map(|candidate| candidate + offset)
        .collect())
}

pub fn candidates_from_regions(path: &Path, contig: &str) -> Result<Vec<i64>> {
    candidates_from_regions_with_offset(path, contig, -1)
}

pub fn candidates_from_regions_cffi(path: &Path, contig: &str) -> Result<Vec<i64>> {
    candidates_from_regions_with_offset(path, contig, -1)
}

#[derive(Debug, Default, Clone)]
struct PosInfo {
    ins_bases: Option<String>,
    ins_length: usize,
    alt_base: i16,
    del_length: usize,
    bq: i8,
    signal_length: i32,
}

#[derive(Debug)]
struct ReadAlignment {
    index: usize,
    record: Record,
    cigar: Vec<bam::record::Cigar>,
    read_start: i64,
    read_end: i64,
    haplotype: usize,
    flanking_start: usize,
    overlap_candidates_num: usize,
    pos_info: Vec<PosInfo>,
}

impl ReadAlignment {
    fn new(record: Record, read_start: i64, read_end: i64) -> Self {
        let cigar = record.cigar().iter().copied().collect();
        Self {
            index: 0,
            record,
            cigar,
            read_start,
            read_end,
            haplotype: 0,
            flanking_start: 0,
            overlap_candidates_num: 0,
            pos_info: vec![],
        }
    }

    fn mq(&self) -> i8 {
        ((self.record.mapq() as i32 * 100 / 60).min(100)) as i8
    }

    fn strand(&self) -> i8 {
        if self.record.flags() & 0x10 == 0 {
            100
        } else {
            50
        }
    }
}

fn compact_record(record: &Record, dwell: bool) -> Result<Record> {
    let cigar = bam::record::CigarString(record.cigar().iter().copied().collect());
    let sequence = record.seq().as_bytes();
    let mut compact = Record::new();
    compact.set(record.qname(), Some(&cigar), &sequence, record.qual());
    compact.set_tid(record.tid());
    compact.set_pos(record.pos());
    compact.set_mapq(record.mapq());
    compact.set_flags(record.flags());
    compact.set_mtid(record.mtid());
    compact.set_mpos(record.mpos());
    compact.set_insert_size(record.insert_size());
    if dwell && let Ok(aux) = record.aux(b"mv") {
        compact.push_aux(b"mv", aux)?;
    }
    Ok(compact)
}

fn haplotype_tag(record: &Record) -> Option<usize> {
    let value = match record.aux(b"HP").ok()? {
        Aux::I8(value) => value as i64,
        Aux::U8(value) => value as i64,
        Aux::I16(value) => value as i64,
        Aux::U16(value) => value as i64,
        Aux::I32(value) => value as i64,
        Aux::U32(value) => value as i64,
        _ => return None,
    };
    Some(value.max(0) as usize)
}

#[derive(Debug, Default)]
struct CandidateAltInfo {
    ins_counter: KhCounter<String>,
    del_counter: KhCounter<usize>,
    acgt_count: [i64; 4],
    depth: i64,
}

#[derive(Debug, Default, Clone)]
struct AltInfoRow {
    alt_base: u8,
    ins_bases: Option<String>,
    del_length: usize,
    has_alt_info: bool,
}

fn acgt_index(base: u8) -> usize {
    match base {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

fn python_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        base @ (b'A' | b'C' | b'G' | b'T') => base,
        _ => b'A',
    }
}

fn python_encoded_base(base: i16) -> u8 {
    if (0..NT16.len() as i16).contains(&base) {
        python_base(NT16[base as usize])
    } else {
        b'A'
    }
}

fn record_alt_count(counts: &mut Vec<(String, i64)>, key: String) {
    if let Some((_, count)) = counts.iter_mut().find(|(current, _)| current == &key) {
        *count += 1;
    } else {
        counts.push((key, 1));
    }
}

fn python_deletion(ref_seq: &[u8], candidate_offset: i64, length: usize) -> String {
    ref_seq
        .get((candidate_offset + 1) as usize..(candidate_offset + 1 + length as i64) as usize)
        .map(String::from_utf8_lossy)
        .unwrap_or_default()
        .to_ascii_uppercase()
}

fn python_alt_counts(
    reads: &[ReadAlignment],
    flanking_indices: &BTreeMap<i64, usize>,
    candidate: i64,
    ref_seq: &[u8],
    reference_start: i64,
    center_ref: u8,
) -> (Vec<(String, i64)>, i64) {
    let mut counts = Vec::new();
    let mut ref_count = 0;
    for read in reads {
        let Some(flanking_index) = flanking_indices.get(&candidate).copied() else {
            continue;
        };
        if flanking_index < read.flanking_start {
            continue;
        }
        let Some(info) = read.pos_info.get(flanking_index - read.flanking_start) else {
            continue;
        };
        if info.alt_base < 0 {
            continue;
        }
        if info.ins_length > 0 {
            if let Some(insertion) = &info.ins_bases {
                record_alt_count(
                    &mut counts,
                    format!(
                        "I{}{}",
                        python_encoded_base(info.alt_base) as char,
                        insertion
                    ),
                );
            }
        } else if info.del_length > 0 {
            let deleted = python_deletion(ref_seq, candidate - reference_start, info.del_length);
            record_alt_count(&mut counts, format!("D{deleted}"));
        } else if info.alt_base >= 0 {
            let base = python_encoded_base(info.alt_base);
            if base == center_ref {
                ref_count += 1;
            } else {
                record_alt_count(&mut counts, format!("X{}", base as char));
            }
        }
    }
    (counts, ref_count)
}

fn position_read<'a>(
    reads: &'a [ReadAlignment],
    names: Option<&BTreeMap<Vec<u8>, Vec<usize>>>,
    read_index: usize,
    flanking_index: usize,
) -> &'a ReadAlignment {
    let mut selected = &reads[read_index];
    let Some(indices) = names.and_then(|names| names.get(selected.record.qname())) else {
        return selected;
    };
    for index in indices {
        let read = &reads[*index];
        if flanking_index >= read.flanking_start
            && flanking_index - read.flanking_start < read.pos_info.len()
        {
            selected = read;
        }
    }
    selected
}

fn channel_base(base: u8) -> i8 {
    match base {
        b'A' | b'a' | b'N' | b'n' => 100,
        b'C' | b'c' => 25,
        b'G' | b'g' => 75,
        b'T' | b't' => 50,
        b'I' | b'i' => -50,
        b'D' | b'd' => -100,
        _ => 0,
    }
}

fn normalize_bq(quality: u8) -> i8 {
    ((quality as i32 * 100 / 40).min(100)) as i8
}

fn normalize_af(value: f32) -> i8 {
    if value < 1.0 {
        (100.0 * value) as i8
    } else {
        100
    }
}

fn read_end(cigar: &[bam::record::Cigar], start: i64) -> i64 {
    start
        + cigar
            .iter()
            .map(|operation| match operation {
                bam::record::Cigar::Match(length)
                | bam::record::Cigar::Equal(length)
                | bam::record::Cigar::Diff(length)
                | bam::record::Cigar::Del(length)
                | bam::record::Cigar::RefSkip(length) => *length as i64,
                _ => 0,
            })
            .sum::<i64>()
}

fn one_deletion_apart(longer: &[u8], shorter: &[u8]) -> bool {
    let prefix = longer
        .iter()
        .zip(shorter)
        .take_while(|(left, right)| left == right)
        .count();
    longer[prefix + 1..] == shorter[prefix..]
}

fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    if a.is_empty() {
        return b.len();
    }
    if b.is_empty() {
        return a.len();
    }
    if a == b {
        return 0;
    }
    if a.len() == b.len()
        && a.iter()
            .zip(b)
            .filter(|(left, right)| left != right)
            .count()
            == 1
    {
        return 1;
    }
    if (a.len() == b.len() + 1 && one_deletion_apart(a, b))
        || (b.len() == a.len() + 1 && one_deletion_apart(b, a))
    {
        return 1;
    }
    let mut cache: Vec<usize> = (1..=a.len()).collect();
    for (b_index, code) in b.iter().enumerate() {
        let mut result = b_index;
        let mut distance = b_index;
        for (index, a_base) in a.iter().enumerate() {
            let b_distance = if code == a_base {
                distance
            } else {
                distance + 1
            };
            distance = cache[index];
            cache[index] = if distance > result {
                if b_distance > result {
                    result + 1
                } else {
                    b_distance
                }
            } else if b_distance > distance {
                distance + 1
            } else {
                b_distance
            };
            result = cache[index];
        }
        let _ = result;
    }
    cache[a.len() - 1]
}

fn slice(sequence: &[u8], start: i64, end: i64) -> Vec<u8> {
    let start = start.clamp(0, sequence.len() as i64) as usize;
    let end = end.clamp(0, sequence.len() as i64) as usize;
    if end <= start {
        vec![]
    } else {
        sequence[start..end].to_vec()
    }
}

fn cigar_prefix_length(
    cigar: &[bam::record::Cigar],
    reference_bases: i64,
    left_cigar_index: usize,
    right_cigar_index: usize,
    consumed: i64,
    reverse: bool,
) -> (i64, i64) {
    let mut ref_pos = 0_i64;
    let mut query_pos = 0_i64;
    for i in left_cigar_index..right_cigar_index {
        let index = if reverse {
            left_cigar_index + right_cigar_index - i - 1
        } else {
            i
        };
        let operation = cigar[index];
        let length = if i == left_cigar_index {
            consumed
        } else {
            operation.len() as i64
        };
        if length == 0 {
            continue;
        }
        match operation {
            bam::record::Cigar::Match(_)
            | bam::record::Cigar::Equal(_)
            | bam::record::Cigar::Diff(_) => {
                query_pos += length;
                ref_pos += length;
                if ref_pos >= reference_bases {
                    return (reference_bases, query_pos + reference_bases - ref_pos);
                }
            }
            bam::record::Cigar::Del(_) => {
                ref_pos += length;
                if ref_pos >= reference_bases {
                    return (reference_bases, query_pos);
                }
            }
            bam::record::Cigar::Ins(_) => query_pos += length,
            bam::record::Cigar::RefSkip(_) => return (reference_bases, query_pos),
            _ => {}
        }
    }
    (0, 0)
}

fn realign_read(
    variant: &PhasedVariant,
    read: &ReadAlignment,
    cigar_index: usize,
    consumed: i64,
    query_pos: i64,
    reference: &[u8],
    ref_start: i64,
) -> i32 {
    let middle_length = read.cigar[cigar_index].len() as i64;
    let left_consumed = consumed.max(0);
    let right_consumed = if consumed < middle_length {
        middle_length - consumed
    } else {
        0
    };
    let (left_ref_bases, left_query_bases) =
        cigar_prefix_length(&read.cigar, 10, 0, cigar_index + 1, left_consumed, true);
    let (right_ref_bases, right_query_bases) = cigar_prefix_length(
        &read.cigar,
        11,
        cigar_index,
        read.cigar.len(),
        right_consumed,
        false,
    );
    let query_start = query_pos - left_query_bases;
    let query_end = query_pos + right_query_bases;
    if query_end == query_start {
        return 0;
    }
    let sequence = read.record.seq();
    let query_start = query_start.clamp(0, sequence.len() as i64) as usize;
    let query_end = query_end.clamp(0, sequence.len() as i64) as usize;
    let query: Vec<u8> = if query_end <= query_start {
        Vec::new()
    } else {
        (query_start..query_end)
            .map(|index| sequence[index])
            .collect()
    };
    let reference_start = variant.position - left_ref_bases - ref_start;
    let reference_end = variant.position + right_ref_bases - ref_start;
    let reference = slice(reference, reference_start, reference_end);
    let mut alt = reference.clone();
    if alt.len() < (left_ref_bases + 1) as usize {
        alt.resize((left_ref_bases + 1) as usize, 0);
    }
    alt[left_ref_bases as usize] = variant.alt_base;
    let distance_ref = levenshtein(&query, &reference);
    let distance_alt = levenshtein(&query, &alt);
    if distance_ref < distance_alt {
        1
    } else if distance_ref > distance_alt {
        2
    } else {
        0
    }
}

fn haplotag_read(
    variants: &[PhasedVariant],
    read: &mut ReadAlignment,
    reference: &[u8],
    ref_start: i64,
) {
    let mut costs: BTreeMap<i32, i32> = BTreeMap::new();
    let mut variant_index = variants.partition_point(|variant| variant.position < read.read_start);
    let mut ref_pos = read.read_start;
    let mut query_pos = 0_i64;
    for (cigar_index, operation) in read.cigar.iter().enumerate() {
        let length = operation.len() as i64;
        match operation {
            bam::record::Cigar::Match(_)
            | bam::record::Cigar::Equal(_)
            | bam::record::Cigar::Diff(_) => {
                while variant_index < variants.len()
                    && variants[variant_index].position < ref_pos + length
                {
                    let variant = &variants[variant_index];
                    let allele = realign_read(
                        variant,
                        read,
                        cigar_index,
                        variant.position - ref_pos,
                        query_pos + variant.position - ref_pos,
                        reference,
                        ref_start,
                    );
                    if allele != 0 {
                        *costs.entry(variant.phase_set).or_default() +=
                            if allele == variant.genotype { 1 } else { -1 };
                    }
                    variant_index += 1;
                }
                query_pos += length;
                ref_pos += length;
            }
            bam::record::Cigar::Ins(_) => {
                if variant_index < variants.len() && variants[variant_index].position == ref_pos {
                    let variant = &variants[variant_index];
                    let allele = realign_read(
                        variant,
                        read,
                        cigar_index,
                        0,
                        query_pos,
                        reference,
                        ref_start,
                    );
                    if allele != 0 {
                        *costs.entry(variant.phase_set).or_default() +=
                            if allele == variant.genotype { 1 } else { -1 };
                    }
                    variant_index += 1;
                }
                query_pos += length;
            }
            bam::record::Cigar::Del(_) => {
                while variant_index < variants.len()
                    && variants[variant_index].position < ref_pos + length
                {
                    let variant = &variants[variant_index];
                    let allele = realign_read(
                        variant,
                        read,
                        cigar_index,
                        variant.position - ref_pos,
                        query_pos,
                        reference,
                        ref_start,
                    );
                    if allele != 0 {
                        *costs.entry(variant.phase_set).or_default() +=
                            if allele == variant.genotype { 1 } else { -1 };
                    }
                    variant_index += 1;
                }
                ref_pos += length;
            }
            bam::record::Cigar::RefSkip(_) => {
                while variant_index < variants.len()
                    && variants[variant_index].position < ref_pos + length
                {
                    variant_index += 1;
                }
                ref_pos += length;
            }
            bam::record::Cigar::SoftClip(_) => query_pos += length,
            _ => {}
        }
    }
    read.read_end = ref_pos;
    if costs.is_empty() || costs.values().all(|value| *value == 0) {
        read.haplotype = 0;
        return;
    }
    let maximum = costs.values().max().copied().unwrap_or(0).max(0);
    let minimum = costs.values().min().copied().unwrap_or(0).min(0);
    read.haplotype = if maximum > minimum.abs() { 1 } else { 2 };
}

fn collect_reads(
    bam_path: &Path,
    reference: &Reference,
    region: &Region,
    candidates: &[i64],
    variants: &[PhasedVariant],
    options: &FullAlignmentOptions,
) -> Result<(Vec<ReadAlignment>, Vec<i64>, Vec<CandidateAltInfo>)> {
    let mut candidate_indices = BTreeMap::new();
    let mut flanking = BTreeSet::new();
    if options.illumina {
        for (index, candidate) in candidates.iter().enumerate() {
            candidate_indices.insert(*candidate, index);
        }
        for candidate in candidates {
            for offset in 0..POSITION_COUNT as i64 {
                flanking.insert((candidate - FLANKING_BASES as i64 + offset).max(0));
            }
        }
    }
    let flanking: Vec<i64> = flanking.into_iter().collect();
    let flanking_indices: BTreeMap<i64, usize> = flanking
        .iter()
        .enumerate()
        .map(|(index, position)| (*position, index))
        .collect();
    let mut candidate_alt: Vec<CandidateAltInfo> = (0..if options.illumina {
        candidates.len()
    } else {
        0
    })
        .map(|_| CandidateAltInfo::default())
        .collect();
    let mut reads = vec![];
    let mut read_names = BTreeSet::new();
    let reference_start = region.start - 2_000_000;
    let ref_seq = reference.fetch_padded(
        &region.contig,
        reference_start,
        region.end + 2_000_000,
        b'A',
    )?;
    let mut reader = IndexedReader::from_path(bam_path)?;
    reader.set_threads(3)?;
    reader.set_reference(reference.path())?;
    reader.fetch((region.contig.as_bytes(), region.start, region.end))?;
    for record in reader.records() {
        let record = record?;
        if record.flags() & FILTER_FLAGS != 0
            || (record.flags() & 0x1 != 0 && record.flags() & 0x2 == 0)
            || record.mapq() < options.min_mq
            || (!options.illumina && !read_names.insert(record.qname().to_vec()))
        {
            continue;
        }
        let start = record.pos();
        let end = read_end(&record.cigar().iter().copied().collect::<Vec<_>>(), start);
        let (flanking_start, overlap) = if options.illumina {
            let start = flanking.partition_point(|position| *position < start);
            let overlap = flanking[start..].partition_point(|position| *position < end);
            (start, overlap)
        } else {
            let next = candidates
                .partition_point(|candidate| *candidate + (FLANKING_BASES as i64) < start);
            if candidates
                .get(next)
                .is_none_or(|candidate| *candidate - FLANKING_BASES as i64 >= end)
            {
                continue;
            }
            (0, 0)
        };
        if options.illumina && overlap == 0 {
            continue;
        }
        let mut read = ReadAlignment::new(record, start, end);
        read.flanking_start = flanking_start;
        read.overlap_candidates_num = overlap;
        read.pos_info = if options.illumina {
            vec![PosInfo::default(); overlap]
        } else {
            Vec::new()
        };
        if !options.force_unphased
            && (options.illumina || options.phasing_info_in_bam)
            && let Some(haplotype) = haplotype_tag(&read.record)
        {
            read.haplotype = haplotype;
        }
        if options.need_haplotagging && read.record.mapq() >= MIN_HAPLOTAG_MQ {
            let tagged = read.haplotype;
            haplotag_read(variants, &mut read, &ref_seq, reference_start);
            read.haplotype = read.haplotype.max(tagged);
        }
        if !options.illumina {
            if options.enable_dwell_time {
                read.record = compact_record(&read.record, true)?;
            }
            read.index = reads.len();
            reads.push(read);
            continue;
        }
        let mut ref_pos = start;
        let mut query_pos = 0_i64;
        let signal_lengths = options
            .enable_dwell_time
            .then(|| signal_lengths(&read.record))
            .flatten();
        for operation in &read.cigar.clone() {
            let length = operation.len() as i64;
            match operation {
                bam::record::Cigar::Match(_)
                | bam::record::Cigar::Equal(_)
                | bam::record::Cigar::Diff(_) => {
                    for position in ref_pos..ref_pos + length {
                        if let Some(flanking_index) = flanking_indices.get(&position).copied()
                            && flanking_index >= flanking_start
                        {
                            let offset = flanking_index - flanking_start;
                            if offset < read.pos_info.len() {
                                let base = read.record.seq().encoded_base(query_pos as usize);
                                read.pos_info[offset].alt_base = base as i16;
                                read.pos_info[offset].bq =
                                    normalize_bq(read.record.qual()[query_pos as usize]);
                                if let Some(signals) = &signal_lengths
                                    && (query_pos as usize) < read.record.seq_len()
                                {
                                    read.pos_info[offset].signal_length =
                                        signals[query_pos as usize];
                                }
                                if let Some(candidate_index) =
                                    candidate_indices.get(&position).copied()
                                {
                                    let base = NT16[base as usize];
                                    let index = acgt_index(base);
                                    candidate_alt[candidate_index].acgt_count[index] += 1;
                                    candidate_alt[candidate_index].depth += 1;
                                }
                            }
                        }
                        query_pos += 1;
                    }
                    ref_pos += length;
                }
                bam::record::Cigar::Del(_) => {
                    if let Some(flanking_index) = flanking_indices.get(&(ref_pos - 1)).copied()
                        && flanking_index >= flanking_start
                    {
                        let offset = flanking_index - flanking_start;
                        if offset < read.pos_info.len() {
                            read.pos_info[offset].del_length = length as usize;
                            if let Some(candidate_index) =
                                candidate_indices.get(&(ref_pos - 1)).copied()
                            {
                                let anchor = read
                                    .record
                                    .seq()
                                    .encoded_base((query_pos - 1).max(0) as usize);
                                let count = &mut candidate_alt[candidate_index].acgt_count
                                    [acgt_index(NT16[anchor as usize])];
                                *count = count.saturating_sub(1);
                                *candidate_alt[candidate_index]
                                    .del_counter
                                    .entry(length as usize) += 1;
                            }
                        }
                    }
                    for position in ref_pos..ref_pos + length {
                        if let Some(flanking_index) = flanking_indices.get(&position).copied()
                            && flanking_index >= flanking_start
                        {
                            let offset = flanking_index - flanking_start;
                            if offset < read.pos_info.len() {
                                read.pos_info[offset].alt_base = -1;
                                if let Some(candidate_index) =
                                    candidate_indices.get(&position).copied()
                                {
                                    candidate_alt[candidate_index].depth += 1;
                                }
                            }
                        }
                    }
                    ref_pos += length;
                }
                bam::record::Cigar::Ins(_) => {
                    if let Some(flanking_index) = flanking_indices.get(&(ref_pos - 1)).copied()
                        && flanking_index >= flanking_start
                    {
                        let offset = flanking_index - flanking_start;
                        if offset < read.pos_info.len() {
                            let mut insertion = String::new();
                            let mut signal_sum = 0;
                            for index in 0..length as usize {
                                let base =
                                    read.record.seq().encoded_base(query_pos as usize + index);
                                insertion.push(NT16[base as usize] as char);
                                if let Some(signals) = &signal_lengths {
                                    signal_sum += signals[query_pos as usize + index];
                                }
                            }
                            read.pos_info[offset].ins_bases = Some(insertion.clone());
                            read.pos_info[offset].ins_length = length as usize;
                            read.pos_info[offset].signal_length += signal_sum;
                            if let Some(candidate_index) =
                                candidate_indices.get(&(ref_pos - 1)).copied()
                            {
                                let anchor = read
                                    .record
                                    .seq()
                                    .encoded_base((query_pos - 1).max(0) as usize);
                                let count = &mut candidate_alt[candidate_index].acgt_count
                                    [acgt_index(NT16[anchor as usize])];
                                *count = count.saturating_sub(1);
                                *candidate_alt[candidate_index].ins_counter.entry(insertion) += 1;
                            }
                        }
                    }
                    query_pos += length;
                }
                bam::record::Cigar::RefSkip(_) => ref_pos += length,
                bam::record::Cigar::SoftClip(_) => query_pos += length,
                _ => {}
            }
        }
        read.index = reads.len();
        reads.push(read);
    }
    Ok((reads, flanking, candidate_alt))
}

struct ReadWindow {
    positions: [PosInfo; POSITION_COUNT],
    covered: [bool; POSITION_COUNT],
}

fn read_window(read: &ReadAlignment, start: i64, dwell: bool) -> ReadWindow {
    let end = start + POSITION_COUNT as i64;
    let mut window = ReadWindow {
        positions: std::array::from_fn(|_| PosInfo::default()),
        covered: [false; POSITION_COUNT],
    };
    let signal_lengths = dwell.then(|| signal_lengths(&read.record)).flatten();
    let mut ref_pos = read.read_start;
    let mut query_pos = 0_i64;
    for operation in &read.cigar {
        if ref_pos > end
            || (ref_pos == end
                && !matches!(
                    operation,
                    bam::record::Cigar::Ins(_) | bam::record::Cigar::Del(_)
                ))
        {
            break;
        }
        let length = operation.len() as i64;
        match operation {
            bam::record::Cigar::Match(_)
            | bam::record::Cigar::Equal(_)
            | bam::record::Cigar::Diff(_) => {
                for position in ref_pos.max(start)..(ref_pos + length).min(end) {
                    let index = (position - start) as usize;
                    let query = (query_pos + position - ref_pos) as usize;
                    let info = &mut window.positions[index];
                    info.alt_base = read.record.seq().encoded_base(query) as i16;
                    info.bq = normalize_bq(read.record.qual()[query]);
                    if let Some(signals) = &signal_lengths {
                        info.signal_length = signals[query];
                    }
                    window.covered[index] = true;
                }
                ref_pos += length;
                query_pos += length;
            }
            bam::record::Cigar::Del(_) => {
                let anchor = ref_pos - 1;
                if anchor >= start && anchor < end && anchor >= read.read_start {
                    window.positions[(anchor - start) as usize].del_length = length as usize;
                }
                for position in ref_pos.max(start)..(ref_pos + length).min(end) {
                    let index = (position - start) as usize;
                    window.positions[index].alt_base = -1;
                    window.covered[index] = true;
                }
                ref_pos += length;
            }
            bam::record::Cigar::Ins(_) => {
                let anchor = ref_pos - 1;
                if anchor >= start && anchor < end && anchor >= read.read_start {
                    let mut insertion = String::new();
                    let mut signal_sum = 0;
                    for index in 0..length as usize {
                        let query = query_pos as usize + index;
                        let base = read.record.seq().encoded_base(query);
                        insertion.push(NT16[base as usize] as char);
                        if let Some(signals) = &signal_lengths {
                            signal_sum += signals[query];
                        }
                    }
                    let info = &mut window.positions[(anchor - start) as usize];
                    info.ins_bases = Some(insertion);
                    info.ins_length = length as usize;
                    info.signal_length += signal_sum;
                }
                query_pos += length;
            }
            bam::record::Cigar::RefSkip(_) => ref_pos += length,
            bam::record::Cigar::SoftClip(_) => query_pos += length,
            _ => {}
        }
    }
    window
}

fn update_candidate_alt(counts: &mut CandidateAltInfo, window: &ReadWindow) {
    if !window.covered[FLANKING_BASES] {
        return;
    }
    let info = &window.positions[FLANKING_BASES];
    if info.alt_base < 0 {
        counts.depth += 1;
    } else {
        let index = acgt_index(NT16[info.alt_base as usize]);
        counts.acgt_count[index] += 1;
        counts.depth += 1;
    }
    if info.del_length > 0 {
        *counts.del_counter.entry(info.del_length) += 1;
    }
    if let Some(insertion) = &info.ins_bases {
        *counts.ins_counter.entry(insertion.clone()) += 1;
    }
}

fn matrix_index(
    candidate: usize,
    depth: usize,
    position: usize,
    channel: usize,
    matrix_depth: usize,
    channels: usize,
) -> usize {
    ((candidate * matrix_depth + depth) * POSITION_COUNT + position) * channels + channel
}

struct PythonRandom {
    state: [u32; 624],
    index: usize,
}

impl PythonRandom {
    fn new(seed: u32) -> Self {
        let mut state = [0_u32; 624];
        state[0] = 19_650_218;
        for index in 1..624 {
            state[index] = 1_812_433_253_u32
                .wrapping_mul(state[index - 1] ^ (state[index - 1] >> 30))
                .wrapping_add(index as u32);
        }
        let mut i = 1_usize;
        let mut j = 0_usize;
        for _ in 0..624 {
            state[i] ^= (state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_mul(1_664_525)
                .wrapping_add(seed)
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= 624 {
                state[0] = state[623];
                i = 1;
            }
            if j >= 1 {
                j = 0;
            }
        }
        for _ in 0..623 {
            state[i] = (state[i]
                ^ (state[i - 1] ^ (state[i - 1] >> 30)).wrapping_mul(1_566_083_941))
            .wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                state[0] = state[623];
                i = 1;
            }
        }
        state[0] = 0x8000_0000;
        Self { state, index: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        if self.index >= 624 {
            for index in 0..624 {
                let value = (self.state[index] & 0x8000_0000)
                    | (self.state[(index + 1) % 624] & 0x7fff_ffff);
                self.state[index] = self.state[(index + 397) % 624]
                    ^ (value >> 1)
                    ^ if value & 1 == 1 { 2_567_483_615 } else { 0 };
            }
            self.index = 0;
        }
        let mut value = self.state[self.index];
        self.index += 1;
        value ^= value >> 11;
        value ^= (value << 7) & 0x9d2c_5680;
        value ^= (value << 15) & 0xefc6_0000;
        value ^= value >> 18;
        value
    }

    fn randbelow(&mut self, n: usize) -> usize {
        let bits = usize::BITS - n.leading_zeros();
        loop {
            let value = (self.next_u32() >> (32 - bits)) as usize;
            if value < n {
                return value;
            }
        }
    }
}

fn sort_reads_python(
    reads: &[ReadAlignment],
    start: i64,
    end: i64,
    matrix_depth: usize,
) -> Vec<isize> {
    let mut read_names = BTreeSet::new();
    let mut overlap: Vec<(usize, usize)> = reads
        [..reads.partition_point(|read| read.read_start < end)]
        .iter()
        .enumerate()
        .filter(|(_, read)| read.read_end > start && read_names.insert(read.record.qname()))
        .map(|(index, read)| (index, read.haplotype))
        .collect();
    let n = overlap.len();
    let read_num = n.min(matrix_depth);
    if n > matrix_depth {
        let setsize = if read_num > 5 {
            21 + 4_usize.pow(((read_num * 3) as f64).log(4.0).ceil() as u32)
        } else {
            21
        };
        let mut random = PythonRandom::new(0);
        let mut sampled = Vec::with_capacity(read_num);
        if n > setsize {
            let mut selected = BTreeSet::new();
            for _ in 0..read_num {
                let mut index = random.randbelow(n);
                while selected.contains(&index) {
                    index = random.randbelow(n);
                }
                selected.insert(index);
                sampled.push(overlap[index]);
            }
        } else {
            for index in 0..read_num {
                let selected = random.randbelow(n - index);
                sampled.push(overlap[selected]);
                overlap[selected] = overlap[n - index - 1];
            }
        }
        sampled.sort_unstable();
        overlap = sampled;
    }
    overlap[..read_num].sort_by_key(|(read_index, haplotype)| (*haplotype, *read_index));
    let padding = matrix_depth - read_num;
    let prefix = padding >> 1;
    let mut indices = vec![-1_isize; matrix_depth];
    for (index, (read_index, _)) in overlap.into_iter().take(read_num).enumerate() {
        indices[prefix + index] = read_index as isize;
    }
    indices
}

fn sort_reads(
    reads: &[ReadAlignment],
    start: i64,
    end: i64,
    matrix_depth: usize,
    illumina: bool,
    random_seed: &mut SamplingRandom,
) -> Vec<isize> {
    if illumina {
        return sort_reads_python(reads, start, end, matrix_depth);
    }
    let mut overlap: Vec<(usize, usize, usize)> = reads
        [..reads.partition_point(|read| read.read_start < end)]
        .iter()
        .enumerate()
        .filter(|(_, read)| read.read_end > start)
        .map(|(index, read)| (index, read.haplotype, read.index))
        .collect();
    let n = overlap.len();
    let read_num = n.min(matrix_depth);
    if n > matrix_depth {
        for i in 0..n - 1 {
            let random = random_seed.next();
            let j = i + random / (libc::RAND_MAX as usize / (n - i) + 1);
            overlap.swap(i, j);
        }
    }
    overlap[..read_num].sort_by_key(|(_, haplotype, global_index)| (*haplotype, *global_index));
    let padding = matrix_depth - read_num;
    let prefix = padding >> 1;
    let mut indices = vec![-1_isize; matrix_depth];
    for (index, (read_index, _, _)) in overlap.into_iter().take(read_num).enumerate() {
        indices[prefix + index] = read_index as isize;
    }
    indices
}

#[allow(clippy::too_many_arguments)]
fn create_full_alignment_impl(
    region: &Region,
    reference_start: i64,
    ref_seq: &[u8],
    candidates: &[i64],
    reads: &[ReadAlignment],
    flanking: Vec<i64>,
    mut candidate_alt: Vec<CandidateAltInfo>,
    options: &FullAlignmentOptions,
    random_seed: &mut SamplingRandom,
    sampled_reads: Option<&[Vec<isize>]>,
) -> Result<FullAlignmentBatch> {
    if candidates.is_empty() {
        return Ok(FullAlignmentBatch {
            shape: vec![
                0,
                options.matrix_depth,
                POSITION_COUNT,
                if options.enable_dwell_time {
                    FULL_ALIGNMENT_CHANNELS + 1
                } else {
                    FULL_ALIGNMENT_CHANNELS
                },
            ],
            ..Default::default()
        });
    }
    let flanking_indices: BTreeMap<i64, usize> = flanking
        .iter()
        .enumerate()
        .map(|(index, position)| (*position, index))
        .collect();
    let read_name_indices = options.illumina.then(|| {
        let mut indices = BTreeMap::<Vec<u8>, Vec<usize>>::new();
        for (index, read) in reads.iter().enumerate() {
            indices
                .entry(read.record.qname().to_vec())
                .or_default()
                .push(index);
        }
        indices
    });
    let channels = if options.enable_dwell_time {
        FULL_ALIGNMENT_CHANNELS + 1
    } else {
        FULL_ALIGNMENT_CHANNELS
    };
    let mut batch = FullAlignmentBatch {
        shape: vec![
            candidates.len(),
            options.matrix_depth,
            POSITION_COUNT,
            channels,
        ],
        values: vec![0; candidates.len() * options.matrix_depth * POSITION_COUNT * channels],
        ..Default::default()
    };

    for (candidate_index, candidate) in candidates.iter().enumerate() {
        let center_ref = ref_seq[(candidate - reference_start) as usize].to_ascii_uppercase();
        let python_alt = options.illumina.then(|| {
            python_alt_counts(
                reads,
                &flanking_indices,
                *candidate,
                ref_seq,
                reference_start,
                center_ref,
            )
        });
        let start = candidate - FLANKING_BASES as i64;
        let end = candidate + FLANKING_BASES as i64 + 1;
        let mut read_windows = Vec::new();
        if !options.illumina {
            let read_end = reads.partition_point(|read| read.read_start < end);
            let mut counts = CandidateAltInfo::default();
            for (read_index, read) in reads[..read_end]
                .iter()
                .enumerate()
                .filter(|(_, read)| read.read_end > start)
            {
                let window = read_window(read, start, options.enable_dwell_time);
                update_candidate_alt(&mut counts, &window);
                read_windows.push((read_index, window));
            }
            candidate_alt.push(counts);
        }
        let read_indices = if let Some(sampled_reads) = sampled_reads {
            sampled_reads[candidate_index].clone()
        } else if *candidate < FLANKING_BASES as i64 {
            vec![-1; options.matrix_depth]
        } else {
            sort_reads(
                reads,
                start,
                end,
                options.matrix_depth,
                options.illumina,
                random_seed,
            )
        };
        let reference_window: Vec<(u8, i8)> = (start..end)
            .map(|position| {
                let base = ref_seq[(position - reference_start) as usize].to_ascii_uppercase();
                (base, channel_base(base))
            })
            .collect();
        let mut alt_rows = vec![AltInfoRow::default(); options.matrix_depth];
        for (depth_index, read_index) in read_indices.into_iter().enumerate() {
            if read_index < 0 {
                continue;
            }
            let haplotype = [60, 30, 90][reads[read_index as usize].haplotype];
            let read_window = (!options.illumina).then(|| {
                let window_index = read_windows
                    .binary_search_by_key(&(read_index as usize), |(index, _)| *index)
                    .unwrap();
                &read_windows[window_index].1
            });
            for (position_index, cp) in (start..end).enumerate() {
                if cp < 0 || (options.illumina && !flanking_indices.contains_key(&cp)) {
                    continue;
                }
                let read = if options.illumina {
                    position_read(
                        reads,
                        read_name_indices.as_ref(),
                        read_index as usize,
                        flanking_indices[&cp],
                    )
                } else {
                    &reads[read_index as usize]
                };
                let info = if options.illumina {
                    let flanking_index = flanking_indices[&cp];
                    let offset = flanking_index as i64 - read.flanking_start as i64;
                    if offset < 0 || offset as usize >= read.overlap_candidates_num {
                        continue;
                    }
                    &read.pos_info[offset as usize]
                } else {
                    &read_window.unwrap().positions[position_index]
                };
                if (!options.illumina && (cp < read.read_start || cp >= read.read_end))
                    || info.alt_base < 0
                {
                    continue;
                }
                let (ref_base, ref_value) = reference_window[position_index];
                let mut alt_value = 0;
                let alt_base = NT16[info.alt_base as usize];
                let is_center = position_index == FLANKING_BASES;
                if info.ins_length > 0 {
                    if let Some(insertion) = &info.ins_bases {
                        let insertion_bytes = insertion.as_bytes();
                        let max_length = insertion_bytes.len().min(POSITION_COUNT - position_index);
                        for (insertion_index, base) in insertion_bytes.iter().enumerate() {
                            if insertion_index < max_length && position_index < POSITION_COUNT - 1 {
                                let index = matrix_index(
                                    candidate_index,
                                    depth_index,
                                    insertion_index + position_index,
                                    6,
                                    options.matrix_depth,
                                    channels,
                                );
                                batch.values[index] = channel_base(*base);
                            }
                        }
                        if is_center {
                            alt_rows[depth_index].alt_base = alt_base;
                            alt_rows[depth_index].ins_bases = Some(insertion.clone());
                            alt_rows[depth_index].has_alt_info = true;
                        }
                    }
                    alt_value = channel_base(b'I');
                } else if info.del_length > 0 {
                    if is_center {
                        alt_rows[depth_index].del_length = info.del_length;
                        alt_rows[depth_index].has_alt_info = true;
                    }
                    alt_value = channel_base(b'D');
                } else if ref_base != alt_base {
                    if is_center {
                        alt_rows[depth_index].alt_base = alt_base;
                        alt_rows[depth_index].has_alt_info = true;
                    }
                    alt_value = if options.illumina
                        && !matches!(
                            alt_base,
                            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
                        ) {
                        0
                    } else {
                        channel_base(alt_base)
                    };
                }
                for (channel, value) in [
                    (0, ref_value),
                    (1, alt_value),
                    (
                        2,
                        if options.illumina && !matches!(alt_base, b'A' | b'C' | b'G' | b'T') {
                            50
                        } else {
                            read.strand()
                        },
                    ),
                    (3, read.mq()),
                    (4, info.bq),
                    (7, haplotype),
                ] {
                    let index = matrix_index(
                        candidate_index,
                        depth_index,
                        position_index,
                        channel,
                        options.matrix_depth,
                        channels,
                    );
                    batch.values[index] = value;
                }
                if options.enable_dwell_time {
                    let index = matrix_index(
                        candidate_index,
                        depth_index,
                        position_index,
                        channels - 1,
                        options.matrix_depth,
                        channels,
                    );
                    batch.values[index] = info.signal_length as i8;
                }
            }
        }

        let candidate_depth = candidate_alt[candidate_index].depth;
        for (depth_index, alt) in alt_rows.iter().enumerate() {
            if !alt.has_alt_info || candidate_depth == 0 {
                continue;
            }
            let af = if let Some((entries, _)) = &python_alt {
                let key = if let Some(insertion) = &alt.ins_bases {
                    format!("I{}{}", python_base(alt.alt_base) as char, insertion)
                } else if alt.del_length > 0 {
                    format!(
                        "D{}",
                        python_deletion(ref_seq, candidate - reference_start, alt.del_length)
                    )
                } else {
                    format!("X{}", python_base(alt.alt_base) as char)
                };
                let count = entries
                    .iter()
                    .find(|(current, _)| current == &key)
                    .map(|(_, count)| *count)
                    .unwrap_or(0);
                normalize_af(count as f32 / candidate_depth as f32)
            } else if let Some(insertion) = &alt.ins_bases {
                let count = candidate_alt[candidate_index]
                    .ins_counter
                    .get(insertion.as_str())
                    .copied()
                    .unwrap_or(0);
                normalize_af(count as f32 / candidate_depth as f32)
            } else if alt.del_length > 0 {
                let count = candidate_alt[candidate_index]
                    .del_counter
                    .get(&alt.del_length)
                    .copied()
                    .unwrap_or(0);
                normalize_af(count as f32 / candidate_depth as f32)
            } else if alt.alt_base != 0 {
                let count = candidate_alt[candidate_index].acgt_count[acgt_index(alt.alt_base)];
                normalize_af(count as f32 / candidate_depth as f32)
            } else {
                0
            };
            if af > 0 {
                for position in 0..POSITION_COUNT {
                    let ref_index = matrix_index(
                        candidate_index,
                        depth_index,
                        position,
                        0,
                        options.matrix_depth,
                        channels,
                    );
                    if batch.values[ref_index] != 0 {
                        let index = matrix_index(
                            candidate_index,
                            depth_index,
                            position,
                            5,
                            options.matrix_depth,
                            channels,
                        );
                        batch.values[index] = af;
                    }
                }
            }
        }

        let raw_ref_count = python_alt.as_ref().map_or_else(
            || candidate_alt[candidate_index].acgt_count[acgt_index(center_ref)],
            |(_, count)| *count,
        );
        let ref_count = if python_alt.is_some() {
            raw_ref_count
        } else {
            let insertions: i64 = candidate_alt[candidate_index]
                .ins_counter
                .iter()
                .map(|(_, count)| *count)
                .sum();
            let deletions: i64 = candidate_alt[candidate_index]
                .del_counter
                .iter()
                .map(|(_, count)| *count)
                .sum();
            raw_ref_count - insertions - deletions
        };
        let mut payload = String::new();
        if let Some((entries, _)) = &python_alt {
            for (allele, count) in entries {
                payload.push_str(&format!("{allele} {count} "));
            }
        } else {
            for (index, base) in b"ACGT".iter().enumerate() {
                if index != acgt_index(center_ref)
                    && candidate_alt[candidate_index].acgt_count[index] > 0
                {
                    payload.push_str(&format!(
                        "X{} {} ",
                        *base as char, candidate_alt[candidate_index].acgt_count[index]
                    ));
                }
            }
            for (insertion, count) in candidate_alt[candidate_index].ins_counter.iter() {
                if options.illumina || insertion.len() <= options.max_indel_length {
                    payload.push_str(&format!("I{}{} {} ", center_ref as char, insertion, count));
                }
            }
            for (length, count) in candidate_alt[candidate_index].del_counter.iter() {
                if options.illumina || *length <= options.max_indel_length {
                    let deleted = ref_seq
                        .get(
                            (candidate - reference_start + 1) as usize
                                ..(candidate - reference_start + 1 + *length as i64) as usize,
                        )
                        .map(String::from_utf8_lossy)
                        .unwrap_or_default();
                    payload.push_str(&format!("D{} {} ", deleted, count));
                }
            }
        }
        if ref_count > 0 {
            payload.push_str(&format!("R{} {} ", center_ref as char, ref_count));
        }
        batch.positions.push(format!(
            "{}:{}:{}",
            region.contig,
            candidate + 1,
            center_ref as char
        ));
        batch
            .alt_info
            .push(format!("{}-{}", candidate_depth, payload.trim_end()));
    }
    Ok(batch)
}

fn overlaps_candidates(candidates: &[i64], start: i64, end: i64) -> bool {
    let next = candidates.partition_point(|candidate| *candidate + (FLANKING_BASES as i64) < start);
    candidates
        .get(next)
        .is_some_and(|candidate| *candidate - (FLANKING_BASES as i64) < end)
}

fn prepare_read(
    record: Record,
    candidates: &[i64],
    variants: &[PhasedVariant],
    options: &FullAlignmentOptions,
    ref_seq: &[u8],
    reference_start: i64,
) -> Result<Option<ReadAlignment>> {
    let start = record.pos();
    let end = read_end(&record.cigar().iter().copied().collect::<Vec<_>>(), start);
    if !overlaps_candidates(candidates, start, end) {
        return Ok(None);
    }
    let mut read = ReadAlignment::new(record, start, end);
    if !options.force_unphased
        && options.phasing_info_in_bam
        && let Some(haplotype) = haplotype_tag(&read.record)
    {
        read.haplotype = haplotype;
    }
    if options.need_haplotagging && read.record.mapq() >= MIN_HAPLOTAG_MQ {
        let tagged = read.haplotype;
        haplotag_read(variants, &mut read, ref_seq, reference_start);
        read.haplotype = read.haplotype.max(tagged);
    }
    if options.enable_dwell_time {
        read.record = compact_record(&read.record, true)?;
    }
    Ok(Some(read))
}

fn full_reference_window(
    reference: &Reference,
    region: &Region,
    options: &FullAlignmentOptions,
) -> Result<(i64, Vec<u8>)> {
    let padding = options
        .max_indel_length
        .saturating_add(1)
        .max(FLANKING_BASES) as i64;
    let start = region.start - padding;
    Ok((
        start,
        reference.fetch_padded(&region.contig, start, region.end + padding, b'A')?,
    ))
}

fn create_full_alignment_chunk(
    reference: &Reference,
    region: &Region,
    candidates: &[i64],
    reads: &[ReadAlignment],
    options: &FullAlignmentOptions,
    random_seed: &mut SamplingRandom,
) -> Result<FullAlignmentBatch> {
    let sampled: Vec<_> = candidates
        .iter()
        .map(|candidate| {
            if *candidate < FLANKING_BASES as i64 {
                vec![-1; options.matrix_depth]
            } else {
                sort_reads(
                    reads,
                    candidate - FLANKING_BASES as i64,
                    candidate + FLANKING_BASES as i64 + 1,
                    options.matrix_depth,
                    false,
                    random_seed,
                )
            }
        })
        .collect();
    let (reference_start, ref_seq) = full_reference_window(reference, region, options)?;
    if candidates.len() < 2 {
        return create_full_alignment_impl(
            region,
            reference_start,
            &ref_seq,
            candidates,
            reads,
            Vec::new(),
            Vec::new(),
            options,
            random_seed,
            Some(&sampled),
        );
    }
    let workers = if candidates.len() >= 4 { 4 } else { 2 };
    let chunk_size = candidates.len().div_ceil(workers);
    let ref_seq = ref_seq.as_slice();
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (candidates, sampled) in candidates
            .chunks(chunk_size)
            .zip(sampled.chunks(chunk_size))
        {
            handles.push(scope.spawn(move || {
                let mut unused_random = SamplingRandom::new(1);
                create_full_alignment_impl(
                    region,
                    reference_start,
                    ref_seq,
                    candidates,
                    reads,
                    Vec::new(),
                    Vec::new(),
                    options,
                    &mut unused_random,
                    Some(sampled),
                )
            }));
        }
        let mut output: Option<FullAlignmentBatch> = None;
        for handle in handles {
            let batch = handle
                .join()
                .map_err(|_| anyhow::anyhow!("full-alignment tensor worker panicked"))??;
            if let Some(output) = &mut output {
                output.append(batch)?;
            } else {
                output = Some(batch);
            }
        }
        output.ok_or_else(|| anyhow::anyhow!("full-alignment tensor workers returned no output"))
    })
}

pub fn create_full_alignment(
    bam_path: &Path,
    reference: &Reference,
    region: &Region,
    candidates: &[i64],
    variants: &[PhasedVariant],
    options: &FullAlignmentOptions,
) -> Result<FullAlignmentBatch> {
    let mut random_seed = SamplingRandom::new(1);
    if options.illumina {
        let (reads, flanking, candidate_alt) =
            collect_reads(bam_path, reference, region, candidates, variants, options)?;
        let (reference_start, ref_seq) = full_reference_window(reference, region, options)?;
        return create_full_alignment_impl(
            region,
            reference_start,
            &ref_seq,
            candidates,
            &reads,
            flanking,
            candidate_alt,
            options,
            &mut random_seed,
            None,
        );
    }
    let channels = if options.enable_dwell_time {
        FULL_ALIGNMENT_CHANNELS + 1
    } else {
        FULL_ALIGNMENT_CHANNELS
    };
    let mut batch = FullAlignmentBatch {
        shape: vec![0, options.matrix_depth, POSITION_COUNT, channels],
        ..Default::default()
    };
    if candidates.is_empty() {
        return Ok(batch);
    }
    let reference_start = region.start - 2_000_000;
    let ref_seq = reference.fetch_padded(
        &region.contig,
        reference_start,
        region.end + 2_000_000,
        b'A',
    )?;
    let mut reader = IndexedReader::from_path(bam_path)?;
    reader.set_threads(3)?;
    reader.set_reference(reference.path())?;
    reader.fetch((region.contig.as_bytes(), region.start, region.end))?;
    let mut records = reader.records();
    let mut read_names = BTreeSet::new();
    let mut next_index = 0;
    let mut pending = Vec::new();
    let mut next_record = None;
    for (chunk_index, chunk) in candidates.chunks(256).enumerate() {
        let start = chunk[0] - FLANKING_BASES as i64;
        let end = chunk[chunk.len() - 1] + FLANKING_BASES as i64 + 1;
        let next_start = candidates
            .get((chunk_index + 1) * 256)
            .map_or(i64::MAX, |candidate| *candidate - FLANKING_BASES as i64);
        pending.retain(|read: &ReadAlignment| read.read_end > start);
        let mut current = std::mem::take(&mut pending);
        let mut buffered = Vec::new();
        loop {
            let next = match next_record.take() {
                Some(record) => Some(Ok(record)),
                None => records.next(),
            };
            let Some(record) = next else {
                break;
            };
            let record = record?;
            if record.pos() >= end {
                next_record = Some(record);
                break;
            }
            if record.flags() & FILTER_FLAGS != 0
                || (record.flags() & 0x1 != 0 && record.flags() & 0x2 == 0)
                || record.mapq() < options.min_mq
                || !read_names.insert(record.qname().to_vec())
            {
                continue;
            }
            buffered.push(record);
        }
        if !buffered.is_empty() {
            let chunk_size = buffered.len().div_ceil(4);
            let mut chunks = Vec::new();
            let mut input = buffered.into_iter();
            for _ in 0..4 {
                let chunk: Vec<_> = input.by_ref().take(chunk_size).collect();
                if !chunk.is_empty() {
                    chunks.push(chunk);
                }
            }
            let ref_seq = ref_seq.as_slice();
            let prepared = std::thread::scope(|scope| {
                let mut handles = Vec::new();
                for records in chunks {
                    handles.push(scope.spawn(move || {
                        records
                            .into_iter()
                            .map(|record| {
                                prepare_read(
                                    record,
                                    candidates,
                                    variants,
                                    options,
                                    ref_seq,
                                    reference_start,
                                )
                            })
                            .collect::<Result<Vec<_>>>()
                    }));
                }
                let mut prepared = Vec::new();
                for handle in handles {
                    prepared.extend(
                        handle
                            .join()
                            .map_err(|_| anyhow::anyhow!("read preparation worker panicked"))??,
                    );
                }
                Ok::<_, anyhow::Error>(prepared)
            })?;
            for mut read in prepared.into_iter().flatten() {
                read.index = next_index;
                next_index += 1;
                if read.read_end > start {
                    current.push(read);
                }
            }
        }
        let chunk_region = Region::new(&region.contig, start.max(0), end)?;
        let next = create_full_alignment_chunk(
            reference,
            &chunk_region,
            chunk,
            &current,
            options,
            &mut random_seed,
        )?;
        batch.append(next)?;
        pending = current
            .into_iter()
            .filter(|read| read.read_end > next_start)
            .collect();
    }
    Ok(batch)
}

pub fn create_candidate_details(
    bam_path: &Path,
    reference: &Reference,
    contig: &str,
    candidates: &[i64],
    illumina: bool,
) -> Result<String> {
    let mut output = String::new();
    for chunk in candidates.chunks(4096) {
        let Some(&first) = chunk.first() else {
            continue;
        };
        let last = *chunk.last().unwrap();
        let contig_length = reference.contig_len(contig) as i64;
        let region = Region::new(
            contig,
            (first - FLANKING_BASES as i64 - 1).max(0),
            (last + FLANKING_BASES as i64 + 2).min(contig_length),
        )?;
        let options = FullAlignmentOptions {
            illumina,
            phasing_info_in_bam: true,
            ..Default::default()
        };
        let (reads, _, _) = collect_reads(bam_path, reference, &region, chunk, &[], &options)?;
        let reference_start = region.start - 2_000_000;
        let reference_sequence = reference.fetch_padded(
            &region.contig,
            reference_start,
            region.end + 2_000_000,
            b'A',
        )?;
        for &candidate in chunk {
            let reference_base =
                reference_sequence[(candidate - reference_start) as usize].to_ascii_uppercase();
            let target = candidate;
            let mut depth = 0_usize;
            let mut alleles = BTreeMap::<String, BTreeSet<(String, usize)>>::new();
            for read in &reads {
                let read_name = String::from_utf8_lossy(read.record.qname()).into_owned();
                let haplotype = read.haplotype;
                let mut allele = None;
                let mut ref_pos = read.read_start;
                let mut query_pos = 0_usize;
                for (index, operation) in read.cigar.iter().enumerate() {
                    let length = operation.len() as i64;
                    match operation {
                        bam::record::Cigar::Match(_)
                        | bam::record::Cigar::Equal(_)
                        | bam::record::Cigar::Diff(_) => {
                            if (ref_pos..ref_pos + length).contains(&target) {
                                let offset = (target - ref_pos) as usize;
                                let base = NT16
                                    [read.record.seq().encoded_base(query_pos + offset) as usize]
                                    .to_ascii_uppercase();
                                let mut indel = None;
                                if ref_pos + length == target + 1 {
                                    match read.cigar.get(index + 1) {
                                        Some(bam::record::Cigar::Ins(size)) => {
                                            let start = query_pos + length as usize;
                                            let sequence: String = (0..*size as usize)
                                                .map(|offset| {
                                                    NT16[read
                                                        .record
                                                        .seq()
                                                        .encoded_base(start + offset)
                                                        as usize]
                                                        as char
                                                })
                                                .collect();
                                            indel = Some(format!("I{}{}", base as char, sequence));
                                        }
                                        Some(bam::record::Cigar::Del(size)) => {
                                            let deleted = reference.fetch(
                                                &region.contig,
                                                ref_pos + length,
                                                ref_pos + length + *size as i64,
                                            )?;
                                            indel = Some(format!(
                                                "D{}",
                                                String::from_utf8_lossy(&deleted)
                                            ));
                                        }
                                        _ => {}
                                    }
                                }
                                allele = Some(if let Some(indel) = indel {
                                    indel
                                } else if base == reference_base {
                                    "R".to_string()
                                } else {
                                    format!("X{}", base as char)
                                });
                            }
                            ref_pos += length;
                            query_pos += length as usize;
                        }
                        bam::record::Cigar::Del(_) => {
                            if (ref_pos..ref_pos + length).contains(&target) {
                                allele = Some("*".to_string());
                            }
                            ref_pos += length;
                        }
                        bam::record::Cigar::RefSkip(_) => ref_pos += length,
                        bam::record::Cigar::SoftClip(_) => query_pos += length as usize,
                        _ => {}
                    }
                }
                if let Some(allele) = allele {
                    depth += 1;
                    alleles
                        .entry(allele)
                        .or_default()
                        .insert((read_name, haplotype));
                }
            }
            let alleles: BTreeMap<String, String> = alleles
                .into_iter()
                .map(|(allele, reads)| {
                    let reads = reads
                        .into_iter()
                        .map(|(name, haplotype)| format!("{name}_{haplotype}"))
                        .collect::<Vec<_>>()
                        .join(" ");
                    (allele, reads)
                })
                .collect();
            output.push_str(&format!(
                "{} {}\t{}\t{}\n",
                region.contig,
                candidate + 1,
                depth,
                serde_json::to_string(&alleles)?
            ));
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_random_matches_cpython() {
        let mut random = PythonRandom::new(0);
        assert_eq!(
            random.state[..10],
            [
                2_147_483_648,
                766_982_754,
                497_961_170,
                3_952_298_588,
                2_331_775_348,
                1_811_986_599,
                3_100_132_149,
                3_188_119_873,
                3_937_547_222,
                215_718_963,
            ]
        );
        let expected = [
            3_626_764_237,
            1_654_615_998,
            3_255_389_356,
            3_823_568_514,
            1_806_341_205,
            173_879_092,
            1_112_038_970,
            4_146_640_122,
            2_195_908_194,
            2_087_043_557,
            1_739_178_872,
            3_943_786_419,
        ];
        for (index, expected) in expected.into_iter().enumerate() {
            assert_eq!(random.next_u32(), expected, "RNG value at index {index}");
        }
    }

    #[test]
    fn python_random_initial_tail_matches_cpython() {
        let random = PythonRandom::new(0);
        assert_eq!(random.state[401], 3_555_715_155);
        assert_eq!(random.state[402], 32_324_691);
    }

    #[test]
    fn python_random_twist_matches_cpython() {
        let mut random = PythonRandom::new(0);
        assert_eq!(random.next_u32(), 3_626_764_237);
        assert_eq!(
            random.state[..12],
            [
                1_372_342_863,
                3_221_959_423,
                4_180_954_279,
                3_990_540_705,
                1_021_773_023,
                2_223_090_966,
                3_843_864_895,
                597_749_037,
                4_038_628_808,
                2_184_695_301,
                2_599_863_370,
                1_285_463_076,
            ]
        );
    }

    #[test]
    fn python_sample_matches_cpython() {
        let mut random = PythonRandom::new(0);
        let population: Vec<usize> = (0..300).collect();
        let mut selected = BTreeSet::new();
        let mut sampled = Vec::new();
        for _ in 0..55 {
            let mut index = random.randbelow(300);
            while selected.contains(&index) {
                index = random.randbelow(300);
            }
            selected.insert(index);
            sampled.push(population[index]);
        }
        assert_eq!(&sampled[..4], &[197, 215, 20, 132]);
        sampled.sort_unstable();
        assert_eq!(
            &sampled[..16],
            &[0, 7, 20, 31, 32, 37, 41, 46, 47, 48, 50, 51, 71, 72, 75, 97]
        );
    }

    #[test]
    fn sampling_stream_matches_platform() {
        let mut random = SamplingRandom::new(1);
        let first = random.next();
        let second = random.next();
        #[cfg(target_os = "linux")]
        assert_eq!((first, second), (1_804_289_383, 846_930_886));
        #[cfg(not(target_os = "linux"))]
        assert_eq!((first, second), (16_807, 282_475_249));
    }

    #[test]
    fn illumina_counts_preserve_anchor_bases_and_skip_deleted_markers() {
        let make_read = |alt_base, insertion: &str| {
            let mut read = ReadAlignment::new(Record::new(), 0, 1);
            read.pos_info = vec![PosInfo {
                ins_bases: Some(insertion.to_string()),
                ins_length: insertion.len(),
                alt_base,
                ..Default::default()
            }];
            read
        };
        let reads = vec![make_read(2, "AA"), make_read(1, "AA"), make_read(-1, "T")];
        let flanking_indices = BTreeMap::from([(0, 0)]);
        let (counts, ref_count) = python_alt_counts(&reads, &flanking_indices, 0, b"CA", 0, b'C');
        assert_eq!(ref_count, 0);
        assert_eq!(
            counts,
            vec![("ICAA".to_string(), 1), ("IAAA".to_string(), 1)]
        );
    }

    #[test]
    fn full_alignment_counts_insertion_after_deletion() {
        let mut record = Record::new();
        record.set(
            b"deletion_insertion",
            Some(&bam::record::CigarString(vec![
                bam::record::Cigar::Match(45),
                bam::record::Cigar::Del(1),
                bam::record::Cigar::Ins(3),
                bam::record::Cigar::Match(45),
            ])),
            &[b'A'; 93],
            &[40; 93],
        );
        record.set_pos(440);
        let read = ReadAlignment::new(record, 440, 531);
        let window = read_window(&read, 469, false);
        let mut counts = CandidateAltInfo::default();
        update_candidate_alt(&mut counts, &window);
        assert_eq!(counts.depth, 1);
        assert_eq!(counts.ins_counter.get("AAA").copied(), Some(1));
    }

    #[test]
    fn read_window_keeps_indel_anchors_at_the_right_edge() {
        let mut insertion = Record::new();
        insertion.set(
            b"insertion",
            Some(&bam::record::CigarString(vec![
                bam::record::Cigar::Match(1),
                bam::record::Cigar::Ins(2),
            ])),
            &[b'A'; 3],
            &[40; 3],
        );
        let read = ReadAlignment::new(insertion, 0, 1);
        let window = read_window(&read, 1 - POSITION_COUNT as i64, false);
        assert_eq!(window.positions[POSITION_COUNT - 1].ins_length, 2);
        assert_eq!(
            window.positions[POSITION_COUNT - 1].ins_bases.as_deref(),
            Some("AA")
        );

        let mut deletion = Record::new();
        deletion.set(
            b"deletion",
            Some(&bam::record::CigarString(vec![
                bam::record::Cigar::Match(1),
                bam::record::Cigar::Del(1),
                bam::record::Cigar::Match(1),
            ])),
            &[b'A'; 3],
            &[40; 3],
        );
        let read = ReadAlignment::new(deletion, 0, 3);
        let window = read_window(&read, 1 - POSITION_COUNT as i64, false);
        assert_eq!(window.positions[POSITION_COUNT - 1].del_length, 1);
    }

    #[test]
    fn realign_read_uses_decoded_sequence_bases() {
        let variant = PhasedVariant {
            position: 110,
            ref_base: b'A',
            alt_base: b'C',
            genotype: 2,
            phase_set: 1,
        };
        let make_read = |sequence: &[u8]| {
            let mut record = Record::new();
            record.set(
                b"read",
                Some(&bam::record::CigarString(vec![bam::record::Cigar::Match(
                    30,
                )])),
                sequence,
                &[40; 30],
            );
            ReadAlignment::new(record, 100, 130)
        };
        let reference = [b'A'; 30];
        let read = make_read(&reference);
        assert_eq!(realign_read(&variant, &read, 0, 10, 10, &reference, 100), 1);
        let mut alternate = reference;
        alternate[10] = b'C';
        let read = make_read(&alternate);
        assert_eq!(realign_read(&variant, &read, 0, 10, 10, &reference, 100), 2);
    }

    #[test]
    fn levenshtein_fast_paths_preserve_edit_distances() {
        for (a, b, expected) in [
            (b"ACGT".as_slice(), b"ACGT".as_slice(), 0),
            (b"ACGT", b"ACTT", 1),
            (b"ACGT", b"AGT", 1),
            (b"ACGT", b"ACGTT", 1),
            (b"ACGT", b"ATGC", 2),
            (b"", b"AC", 2),
        ] {
            assert_eq!(levenshtein(a, b), expected);
            assert_eq!(levenshtein(b, a), expected);
        }
    }

    #[test]
    fn candidate_details_use_the_candidate_base_and_indel_anchor() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data");
        let reference = Reference::open(&data.join("toy.fa")).unwrap();
        let details = create_candidate_details(
            &data.join("toy.bam"),
            &reference,
            "chrToy",
            &[499, 699, 899, 1099],
            false,
        )
        .unwrap();
        for (position, allele, expected_reads) in [
            (500, "XC", 6),
            (500, "R", 6),
            (700, "ITGGG", 12),
            (900, "DCA", 12),
            (1100, "XC", 4),
            (1100, "XG", 4),
            (1100, "R", 4),
        ] {
            let row = details
                .lines()
                .find(|line| line.starts_with(&format!("chrToy {position}\t")))
                .unwrap();
            let fields: Vec<&str> = row.split('\t').collect();
            assert_eq!(fields[1], "12", "depth at {position}");
            let reads: BTreeMap<String, String> = serde_json::from_str(fields[2]).unwrap();
            assert_eq!(
                reads
                    .get(allele)
                    .map(|names| names.split_whitespace().count()),
                Some(expected_reads),
                "allele {allele} at {position}"
            );
        }
    }
}
