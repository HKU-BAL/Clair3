use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

use super::{kh_counter::KhCounter, raw_pileup};
use anyhow::{Result, bail};
use rust_htslib::bam::{self, IndexedReader, pileup::Indel};

use crate::{
    config::{FLANKING_BASES, PILEUP_CHANNELS, POSITION_COUNT},
    io::{
        interval::{IntervalMap, Region},
        reference::Reference,
    },
};

const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const FILTER_FLAGS: u16 = 0x4 | 0x8 | 0x100 | 0x800;

#[derive(Debug, Clone)]
pub struct PileupOptions<'a> {
    pub min_depth: i64,
    pub min_snp_af: f64,
    pub min_indel_af: f64,
    pub min_mq: u8,
    pub max_indel_length: usize,
    pub call_snp_only: bool,
    pub call_head_tail: bool,
    pub exclude_deletions_from_af_denominator: bool,
    pub allow_empty_flanking: bool,
    pub fast_mode: bool,
    pub gvcf: bool,
    pub illumina: bool,
    pub confident_regions: Option<&'a IntervalMap>,
    pub known_positions: Option<&'a BTreeSet<i64>>,
}

impl Default for PileupOptions<'_> {
    fn default() -> Self {
        Self {
            min_depth: 2,
            min_snp_af: 0.08,
            min_indel_af: 0.15,
            min_mq: 5,
            max_indel_length: 50,
            call_snp_only: false,
            call_head_tail: false,
            exclude_deletions_from_af_denominator: false,
            allow_empty_flanking: false,
            fast_mode: false,
            gvcf: false,
            illumina: false,
            confident_regions: None,
            known_positions: None,
        }
    }
}

#[derive(Debug, Default)]
pub struct PileupBatch {
    pub shape: Vec<usize>,
    pub values: Vec<i32>,
    pub positions: Vec<String>,
    pub alt_info: Vec<String>,
    pub gvcf_ref_counts: Vec<u64>,
    pub gvcf_total_counts: Vec<u64>,
    pub has_pileup_rows: bool,
}

impl PileupBatch {
    pub fn values_i8(&self) -> Vec<i8> {
        self.values.iter().map(|value| *value as i8).collect()
    }

    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

#[derive(Debug)]
struct Candidate {
    position: i64,
    ref_base: u8,
    alt_payload: String,
}

#[derive(Debug)]
struct PositionCounts {
    position: i64,
    matrix: [i32; PILEUP_CHANNELS],
    candidate: Option<Candidate>,
    gvcf_ref: u64,
    gvcf_total: u64,
}

fn base_index(encoded: u8, reverse: bool) -> Option<usize> {
    let code = encoded & 0x0f;
    let index = match code {
        1 => 0,
        2 => 1,
        4 => 2,
        8 => 3,
        _ => return None,
    };
    Some(if reverse { index + 9 } else { index })
}

fn reference_index(base: u8) -> usize {
    match base {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

fn insertion_stats(counts: &BTreeMap<String, i64>) -> (i64, i64) {
    counts.values().fold((0, 0), |(sum, maximum), value| {
        (sum + *value, maximum.max(*value))
    })
}

fn record_alt(alt_order: &mut Vec<String>, key: String) {
    if !alt_order.contains(&key) {
        alt_order.push(key);
    }
}

#[allow(clippy::too_many_arguments)]
fn passes_non_ref_majority(
    ref_base: u8,
    major_alt_base: u8,
    ref_count: i64,
    alt_count: i64,
    del_count: i64,
    ins_count: i64,
    illumina: bool,
    pileup_first_is_ref: bool,
) -> bool {
    if illumina {
        return !pileup_first_is_ref;
    }
    ref_count < alt_count
        || ref_count < del_count
        || ref_count < ins_count
        || (ref_count > 0 && ref_count == alt_count && ref_base < major_alt_base)
}

fn deletion_parts(key: &str) -> (usize, &str) {
    key[1..].split_once(':').map_or((0, ""), |(length, bases)| {
        (length.parse().unwrap_or(0), bases)
    })
}

#[allow(clippy::too_many_arguments)]
fn alternate_payload(
    ref_base: u8,
    ref_count: i64,
    matrix: &[i32; PILEUP_CHANNELS],
    deletions_f: &BTreeMap<usize, i64>,
    deletions_r: &BTreeMap<usize, i64>,
    insertions: &KhCounter<String>,
    alt_order: &[String],
    max_indel_length: usize,
    illumina: bool,
) -> String {
    let mut output = String::new();
    let mut ref_depth = ref_count;
    if illumina {
        for allele in alt_order {
            match allele.as_bytes()[0] {
                b'X' => {
                    let base = allele.as_bytes()[1];
                    let index = reference_index(base);
                    let count = matrix[index] + matrix[index + 9];
                    if count > 0 {
                        output.push_str(&format!("X{} {} ", base as char, count));
                    }
                }
                b'D' => {
                    let (length, bases) = deletion_parts(allele);
                    let count = deletions_f.get(&length).copied().unwrap_or(0)
                        + deletions_r.get(&length).copied().unwrap_or(0);
                    ref_depth -= count;
                    if count > 0 && length <= max_indel_length {
                        output.push_str(&format!("D{bases} {count} "));
                    }
                }
                b'I' => {
                    let insertion = &allele[2..];
                    let count = insertions.get(insertion).copied().unwrap_or(0);
                    ref_depth -= count;
                    if insertion.len() <= max_indel_length {
                        output.push_str(&format!("{allele} {count} "));
                    }
                }
                _ => {}
            }
        }
        if ref_depth > 0 {
            output.push_str(&format!("R{} {} ", ref_base as char, ref_depth));
        }
        return output.trim_end().to_string();
    }
    for &base in b"ACGT" {
        if base == ref_base {
            continue;
        }
        let index = reference_index(base);
        let count = matrix[index] + matrix[index + 9];
        if count > 0 {
            output.push_str(&format!("X{} {} ", base as char, count));
        }
    }
    let mut deletions: Vec<_> = alt_order
        .iter()
        .filter(|key| key.starts_with('D'))
        .collect();
    deletions.sort_by_key(|key| deletion_parts(key).0);
    for key in deletions {
        let (length, bases) = deletion_parts(key);
        let count = deletions_f.get(&length).copied().unwrap_or(0)
            + deletions_r.get(&length).copied().unwrap_or(0);
        ref_depth -= count;
        if count > 0 && length <= max_indel_length {
            output.push_str(&format!("D{bases} {count} "));
        }
    }
    for (insertion, count) in insertions.iter() {
        ref_depth -= count;
        if insertion.len() <= max_indel_length {
            output.push_str(&format!("I{}{} {} ", ref_base as char, insertion, count));
        }
    }
    if ref_depth > 0 {
        output.push_str(&format!("R{} {} ", ref_base as char, ref_depth));
    }
    output.trim_end().to_string()
}

fn deleted_anchor_query_offset(record: &raw_pileup::RawRecord, anchor: i64) -> Option<usize> {
    let mut reference_position = record.pos();
    let mut query_position = 0_usize;
    for (index, operation) in record.cigar().iter().enumerate() {
        let length = (operation >> 4) as i64;
        match operation & 0x0f {
            2 => {
                if reference_position <= anchor
                    && anchor < reference_position + length
                    && record
                        .cigar()
                        .get(index + 1)
                        .is_some_and(|operation| operation & 0x0f == 1)
                {
                    return Some(query_position);
                }
                reference_position += length;
            }
            0 | 7 | 8 => {
                reference_position += length;
                query_position += length as usize;
            }
            1 => query_position += length as usize,
            3 => reference_position += length,
            4 => query_position += length as usize,
            _ => {}
        }
    }
    None
}

#[allow(clippy::map_entry)]
fn process_position(
    pileup: &raw_pileup::Pileup,
    region: &Region,
    reference: &Reference,
    ref_seq: &[u8],
    options: &PileupOptions,
    previous_position: &mut Option<i64>,
    contiguous: &mut usize,
) -> Result<Option<PositionCounts>> {
    let position = pileup.position as i64;
    if position < region.start || position >= region.end {
        return Ok(None);
    }
    let current_contiguous = match *previous_position {
        Some(previous) if previous + 1 == position => *contiguous + 1,
        _ => 0,
    };

    let mut matrix = [0_i32; PILEUP_CHANNELS];
    let mut deletions_f: BTreeMap<usize, i64> = BTreeMap::new();
    let mut deletions_r: BTreeMap<usize, i64> = BTreeMap::new();
    let mut insertions_f: BTreeMap<String, i64> = BTreeMap::new();
    let mut insertions_r: BTreeMap<String, i64> = BTreeMap::new();
    let mut insertions_all: KhCounter<String> = KhCounter::default();
    let ref_offset = (position - region.start) as usize;
    let ref_base = ref_seq
        .get(ref_offset)
        .copied()
        .unwrap_or(b'A')
        .to_ascii_uppercase();
    let ref_index = reference_index(ref_base);
    let mut alt_order = Vec::new();
    let mut first_seen = BTreeMap::new();
    let mut next_order = 0_usize;
    let mut depth = 0_i64;
    let mut deleted_bases = 0_i64;
    let mut any_valid = false;

    for raw_alignment in pileup.alignments {
        let alignment = bam::pileup::Alignment::new(raw_alignment);
        let record = raw_pileup::RawRecord::new(raw_alignment.b);
        if record.mapq() < options.min_mq
            || record.flags() & FILTER_FLAGS != 0
            || (record.flags() & 0x1 != 0 && record.flags() & 0x2 == 0)
        {
            continue;
        }
        any_valid = true;
        if alignment.is_refskip() {
            continue;
        }
        let reverse = record.flags() & 0x10 != 0;
        let is_del = alignment.is_del();
        if is_del {
            depth += 1;
            deleted_bases += 1;
            matrix[if reverse { 17 } else { 8 }] += 1;
        } else if let Some(qpos) = alignment.qpos() {
            let encoded = record.encoded_base(qpos);
            if let Some(index) = base_index(encoded, reverse) {
                depth += 1;
                matrix[index] += 1;
                let base = NT16[encoded as usize];
                if !first_seen.contains_key(&base) {
                    first_seen.insert(base, next_order);
                    next_order += 1;
                }
                if base != ref_base {
                    record_alt(&mut alt_order, format!("X{}", base as char));
                }
            }
        }

        match alignment.indel() {
            Indel::Del(length) => {
                if length > 0 {
                    let target = if reverse {
                        &mut deletions_r
                    } else {
                        &mut deletions_f
                    };
                    *target.entry(length as usize).or_default() += 1;
                    let deleted = if ref_offset + 1 + length as usize <= ref_seq.len() {
                        String::from_utf8_lossy(
                            &ref_seq[ref_offset + 1..ref_offset + 1 + length as usize],
                        )
                        .into_owned()
                    } else {
                        String::from_utf8(reference.fetch_padded(
                            &region.contig,
                            position + 1,
                            position + 1 + length as i64,
                            b'A',
                        )?)?
                    };
                    if !first_seen.contains_key(&b'D') {
                        first_seen.insert(b'D', next_order);
                        next_order += 1;
                    }
                    record_alt(&mut alt_order, format!("D{length}:{deleted}"));
                }
            }
            Indel::Ins(length) => {
                if !first_seen.contains_key(&b'I') {
                    first_seen.insert(b'I', next_order);
                    next_order += 1;
                }
                let query_position = if is_del {
                    deleted_anchor_query_offset(&record, pileup.position as i64)
                } else {
                    alignment.qpos()
                };
                if let Some(query_position) = query_position {
                    let first = if is_del { 0 } else { 1 };
                    let mut insertion = String::new();
                    for index in first..length as usize + first {
                        let source = query_position + index;
                        if source >= record.seq_len() {
                            bail!("insertion exceeds read length at {}", pileup.position)
                        }
                        insertion.push(NT16[record.encoded_base(source) as usize] as char);
                    }
                    if reverse {
                        *insertions_r.entry(insertion.clone()).or_default() += 1;
                    } else {
                        *insertions_f.entry(insertion.clone()).or_default() += 1;
                    }
                    *insertions_all.entry(insertion.clone()) += 1;
                    record_alt(
                        &mut alt_order,
                        format!("I{}{}", ref_base as char, insertion),
                    );
                }
            }
            Indel::None => {}
        }
    }
    if !any_valid {
        return Ok(None);
    }
    *contiguous = current_contiguous;
    *previous_position = Some(position);

    for (length, count) in &deletions_f {
        matrix[6] += *count as i32;
        let _ = length;
    }
    if let Some(best) = deletions_f.values().max() {
        matrix[7] = *best as i32;
    }
    for (length, count) in &deletions_r {
        matrix[15] += *count as i32;
        let _ = length;
    }
    if let Some(best) = deletions_r.values().max() {
        matrix[16] = *best as i32;
    }
    let (ins_f_sum, ins_f_best) = insertion_stats(&insertions_f);
    let (ins_r_sum, ins_r_best) = insertion_stats(&insertions_r);
    matrix[4] = ins_f_sum as i32;
    matrix[5] = ins_f_best as i32;
    matrix[13] = ins_r_sum as i32;
    matrix[14] = ins_r_best as i32;

    let forward_sum: i64 = (0..4).map(|index| matrix[index] as i64).sum();
    let reverse_sum: i64 = (0..4).map(|index| matrix[index + 9] as i64).sum();
    let ref_count = matrix[ref_index] as i64 + matrix[ref_index + 9] as i64;
    let mut alt_count = 0_i64;
    let mut major_alt_base = 0_u8;
    let mut all_alt_count = 0_u64;
    let del_count: i64 = deletions_f.values().sum::<i64>() + deletions_r.values().sum::<i64>();
    let ins_count = ins_f_sum + ins_r_sum;
    let mut pileup_list = Vec::new();
    for (index, base) in b"ACGT".iter().enumerate() {
        let current = matrix[index] as i64 + matrix[index + 9] as i64;
        if current > 0 {
            pileup_list.push((*base, current));
        }
        if index != ref_index && current > alt_count {
            major_alt_base = *base;
            all_alt_count += current as u64;
            alt_count = current;
        }
    }
    if ins_count > 0 {
        pileup_list.push((b'I', ins_count));
    }
    if del_count > 0 {
        pileup_list.push((b'D', del_count));
    }
    pileup_list.sort_by(|left, right| {
        right
            .1
            .cmp(&left.1)
            .then_with(|| first_seen.get(&left.0).cmp(&first_seen.get(&right.0)))
    });
    let first_is_ref = pileup_list
        .first()
        .is_none_or(|(base, _)| *base == ref_base);
    let non_ref_majority = passes_non_ref_majority(
        ref_base,
        major_alt_base,
        ref_count,
        alt_count,
        del_count,
        ins_count,
        options.illumina,
        first_is_ref,
    );
    matrix[ref_index] = -(forward_sum as i32);
    matrix[ref_index + 9] = -(reverse_sum as i32);

    let effective_depth = if options.exclude_deletions_from_af_denominator {
        (depth - deleted_bases).max(1)
    } else {
        depth.max(1)
    };
    let snp_af = if options.fast_mode {
        options.min_snp_af.max(0.15)
    } else {
        options.min_snp_af
    };
    let pass_snp_af = alt_count as f64 / effective_depth as f64 >= snp_af
        && (!options.fast_mode || alt_count >= 4);
    let mut pass_af = if options.call_snp_only {
        pass_snp_af
    } else {
        non_ref_majority
            || pass_snp_af
            || del_count as f64 / effective_depth as f64 >= options.min_indel_af
            || ins_count as f64 / effective_depth as f64 >= options.min_indel_af
    };
    if depth > 0 && (options.min_snp_af == 0.0 || options.min_indel_af == 0.0) {
        pass_af = true;
    }
    pass_af &=
        effective_depth >= options.min_depth && matches!(ref_base, b'A' | b'C' | b'G' | b'T');
    if !options.call_head_tail && !options.allow_empty_flanking {
        pass_af &= *contiguous >= FLANKING_BASES;
    }

    let candidate = pass_af.then(|| {
        let payload = alternate_payload(
            ref_base,
            ref_count,
            &matrix,
            &deletions_f,
            &deletions_r,
            &insertions_all,
            &alt_order,
            options.max_indel_length,
            options.illumina,
        );
        Candidate {
            position,
            ref_base,
            alt_payload: format!("{}-{}", depth, payload),
        }
    });
    let gvcf_total = ref_count as u64 + all_alt_count + del_count as u64 + ins_count as u64;
    Ok(Some(PositionCounts {
        position,
        matrix,
        candidate,
        gvcf_ref: ref_count as u64,
        gvcf_total,
    }))
}

fn append_candidate(
    batch: &mut PileupBatch,
    contig: &str,
    candidate: &Candidate,
    tensor: &[[i32; PILEUP_CHANNELS]],
) {
    batch
        .values
        .extend(tensor.iter().flat_map(|row| row.iter().copied()));
    batch.positions.push(format!(
        "{}:{}:{}",
        contig,
        candidate.position + 1,
        candidate.ref_base as char
    ));
    batch.alt_info.push(candidate.alt_payload.clone());
}

fn finalize_candidates(
    batch: &mut PileupBatch,
    contig: &str,
    rows: &[PositionCounts],
    options: &PileupOptions,
    pad_head: bool,
    pad_tail: bool,
) {
    let mut chunks: Vec<&[PositionCounts]> = vec![];
    let mut start = 0;
    for index in 1..=rows.len() {
        if index == rows.len() || rows[index].position != rows[index - 1].position + 1 {
            chunks.push(&rows[start..index]);
            start = index;
        }
    }
    for row in rows {
        let Some(candidate) = &row.candidate else {
            continue;
        };
        let position = candidate.position + 1;
        if options
            .known_positions
            .is_some_and(|positions| !positions.contains(&position))
        {
            continue;
        }
        if options
            .confident_regions
            .is_some_and(|tree| !tree.overlaps(contig, position - 1, position + 1))
        {
            continue;
        }
        for chunk in &chunks {
            let first = chunk[0].position;
            let last = chunk[chunk.len() - 1].position;
            let window_start = position - FLANKING_BASES as i64;
            let window_end = position + FLANKING_BASES as i64 + 1;
            let offset = window_start - first - 1;
            if window_start > first
                && window_end <= last
                && offset >= 0
                && offset + POSITION_COUNT as i64 <= chunk.len() as i64
            {
                let source = &chunk[offset as usize..offset as usize + POSITION_COUNT];
                if !options.allow_empty_flanking
                    && source
                        .iter()
                        .any(|row| row.matrix.iter().all(|value| *value == 0))
                {
                    continue;
                }
                let tensor: Vec<[i32; PILEUP_CHANNELS]> =
                    source.iter().map(|row| row.matrix).collect();
                append_candidate(batch, contig, candidate, &tensor);
            }
            if position < first || position > last {
                continue;
            }
            if pad_head && window_start <= first {
                let take = offset + POSITION_COUNT as i64;
                if take < 0 || take as usize > chunk.len() {
                    continue;
                }
                let mut tensor = vec![[0_i32; PILEUP_CHANNELS]; (-offset) as usize];
                tensor.extend(chunk[..take as usize].iter().map(|row| row.matrix));
                if tensor.len() == POSITION_COUNT {
                    append_candidate(batch, contig, candidate, &tensor);
                }
            }
            if pad_tail && window_end > last && offset >= 0 && offset as usize <= chunk.len() {
                let padding = window_end - last - 2;
                if padding > 0 {
                    let mut tensor: Vec<[i32; PILEUP_CHANNELS]> = chunk[offset as usize..]
                        .iter()
                        .map(|row| row.matrix)
                        .collect();
                    tensor.extend(std::iter::repeat_n(
                        [0_i32; PILEUP_CHANNELS],
                        padding as usize,
                    ));
                    if tensor.len() == POSITION_COUNT {
                        append_candidate(batch, contig, candidate, &tensor);
                    }
                } else if offset + POSITION_COUNT as i64 <= chunk.len() as i64 {
                    let tensor: Vec<[i32; PILEUP_CHANNELS]> = chunk
                        [offset as usize..offset as usize + POSITION_COUNT]
                        .iter()
                        .map(|row| row.matrix)
                        .collect();
                    append_candidate(batch, contig, candidate, &tensor);
                }
            }
        }
    }
}

fn create_pileup_impl(
    bam_path: &Path,
    reference: &Reference,
    region: &Region,
    options: &PileupOptions,
    pad_head: bool,
    pad_tail: bool,
) -> Result<PileupBatch> {
    let ref_seq = reference.fetch_region(region)?;
    let mut reader = IndexedReader::from_path(bam_path)?;
    reader.set_reference(reference.path())?;
    reader.fetch((region.contig.as_bytes(), region.start, region.end))?;
    let mut pileups = raw_pileup::Pileups::new(&mut reader)?;
    pileups.set_max_depth(1 << 20);
    let mut rows = vec![];
    let mut previous_position = None;
    let mut contiguous = 0;
    for pileup in pileups {
        let pileup = pileup?;
        if let Some(row) = process_position(
            &pileup,
            region,
            reference,
            &ref_seq,
            options,
            &mut previous_position,
            &mut contiguous,
        )? {
            rows.push(row);
        }
    }

    let mut batch = PileupBatch {
        gvcf_ref_counts: vec![0; (region.end - region.start) as usize],
        gvcf_total_counts: vec![0; (region.end - region.start) as usize],
        has_pileup_rows: !rows.is_empty(),
        ..Default::default()
    };
    if options.gvcf {
        for row in &rows {
            let index = (row.position - region.start) as usize;
            batch.gvcf_ref_counts[index] = row.gvcf_ref;
            batch.gvcf_total_counts[index] = row.gvcf_total;
        }
    } else {
        batch.gvcf_ref_counts.clear();
        batch.gvcf_total_counts.clear();
    }
    finalize_candidates(
        &mut batch,
        &region.contig,
        &rows,
        options,
        pad_head,
        pad_tail,
    );
    batch.shape = vec![batch.positions.len(), POSITION_COUNT, PILEUP_CHANNELS];
    Ok(batch)
}

fn append_pileup_batch(
    output: &mut PileupBatch,
    chunk: PileupBatch,
    expanded_start: i64,
    central_start: i64,
    central_end: i64,
) -> Result<()> {
    let sample_size = POSITION_COUNT * PILEUP_CHANNELS;
    for (index, position) in chunk.positions.iter().enumerate() {
        let position_number: i64 = position.rsplit(':').nth(1).unwrap_or("").parse()?;
        if position_number <= central_start || position_number > central_end {
            continue;
        }
        output
            .values
            .extend_from_slice(&chunk.values[index * sample_size..(index + 1) * sample_size]);
        output.positions.push(position.clone());
        output.alt_info.push(chunk.alt_info[index].clone());
    }
    output.has_pileup_rows |= chunk.has_pileup_rows;
    if !chunk.gvcf_ref_counts.is_empty() {
        let start = (central_start - expanded_start) as usize;
        let end = (central_end - expanded_start) as usize;
        output
            .gvcf_ref_counts
            .extend_from_slice(&chunk.gvcf_ref_counts[start..end]);
        output
            .gvcf_total_counts
            .extend_from_slice(&chunk.gvcf_total_counts[start..end]);
    }
    Ok(())
}

fn create_pileup_chunked(
    bam_path: &Path,
    reference: &Reference,
    region: &Region,
    options: &PileupOptions,
    chunk_size: i64,
) -> Result<PileupBatch> {
    if chunk_size < POSITION_COUNT as i64 {
        bail!("pileup chunk size must be at least {POSITION_COUNT}");
    }
    if region.start == region.end {
        return create_pileup_impl(
            bam_path,
            reference,
            region,
            options,
            options.call_head_tail,
            options.call_head_tail,
        );
    }
    let contig_length = reference.contig_len(&region.contig) as i64;
    let mut output = PileupBatch {
        shape: vec![0, POSITION_COUNT, PILEUP_CHANNELS],
        ..Default::default()
    };
    let mut start = region.start;
    while start < region.end {
        let end = (start + chunk_size).min(region.end);
        let expanded = Region::new(
            &region.contig,
            (start - 33).max(0),
            (end + 33).min(contig_length),
        )?;
        let chunk = create_pileup_impl(
            bam_path,
            reference,
            &expanded,
            options,
            options.call_head_tail && start == region.start,
            options.call_head_tail && end == region.end,
        )?;
        append_pileup_batch(&mut output, chunk, expanded.start, start, end)?;
        start = end;
    }
    output.shape[0] = output.positions.len();
    Ok(output)
}

pub fn create_pileup(
    bam_path: &Path,
    reference: &Reference,
    region: &Region,
    options: &PileupOptions,
) -> Result<PileupBatch> {
    create_pileup_chunked(bam_path, reference, region, options, 262_144)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunked_pileup_matches_whole_region() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let data = root.join("tests/fixtures/data");
        let reference = Reference::open(&data.join("toy.fa")).unwrap();
        let region = Region::new("chrToy", 0, 4096).unwrap();
        let options = PileupOptions {
            call_head_tail: true,
            gvcf: true,
            ..Default::default()
        };
        let expected = create_pileup_impl(
            &data.join("toy.bam"),
            &reference,
            &region,
            &options,
            true,
            true,
        )
        .unwrap();
        for chunk_size in [34, 64, 127] {
            let actual = create_pileup_chunked(
                &data.join("toy.bam"),
                &reference,
                &region,
                &options,
                chunk_size,
            )
            .unwrap();
            assert_eq!(actual.shape, expected.shape, "shape at {chunk_size}");
            assert_eq!(actual.values, expected.values, "values at {chunk_size}");
            assert_eq!(
                actual.positions, expected.positions,
                "positions at {chunk_size}"
            );
            assert_eq!(
                actual.alt_info, expected.alt_info,
                "alternate information at {chunk_size}"
            );
            assert_eq!(actual.gvcf_ref_counts, expected.gvcf_ref_counts);
            assert_eq!(actual.gvcf_total_counts, expected.gvcf_total_counts);
            assert_eq!(actual.has_pileup_rows, expected.has_pileup_rows);
        }
    }

    #[test]
    fn native_tie_break_and_illumina_first_observed_rules_differ() {
        let tied = |reference, alternate, illumina, first_is_ref| {
            passes_non_ref_majority(reference, alternate, 1, 1, 0, 0, illumina, first_is_ref)
        };
        assert!(tied(b'A', b'G', false, true));
        assert!(!tied(b'T', b'A', false, false));
        assert!(tied(b'T', b'A', true, false));
        assert!(!tied(b'A', b'G', true, true));
    }
}
