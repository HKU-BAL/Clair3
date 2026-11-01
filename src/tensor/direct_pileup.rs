use std::collections::BTreeMap;

use anyhow::{Result, bail};
use rust_htslib::bam::{IndexedReader, Read as BamRead, Record, record::Cigar};

use super::{
    Candidate, FLANKING_BASES, PILEUP_CHANNELS, POSITION_COUNT, PileupBatch, PileupOptions,
    PositionCounts, add_deletion, add_insertion, alternate_payload, base_index,
    finalize_candidates, insertion_stats, passes_non_ref_majority,
};
use crate::{
    io::{interval::Region, reference::Reference},
    tensor::kh_counter::KhCounter,
};

const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

struct RawPosition {
    matrix: [i32; PILEUP_CHANNELS],
    deletions_f: Vec<(usize, i64)>,
    deletions_r: Vec<(usize, i64)>,
    insertions_f: Vec<(String, i64)>,
    insertions_r: Vec<(String, i64)>,
    insertions_all: KhCounter<String>,
    alt_order: Vec<String>,
    first_seen: [u32; 18],
    next_order: u32,
    depth: i64,
    deleted_bases: i64,
    forward_depth: i64,
    reverse_depth: i64,
    any_valid: bool,
}

impl Default for RawPosition {
    fn default() -> Self {
        Self {
            matrix: [0; PILEUP_CHANNELS],
            deletions_f: Vec::new(),
            deletions_r: Vec::new(),
            insertions_f: Vec::new(),
            insertions_r: Vec::new(),
            insertions_all: KhCounter::default(),
            alt_order: Vec::new(),
            first_seen: [u32::MAX; 18],
            next_order: 0,
            depth: 0,
            deleted_bases: 0,
            forward_depth: 0,
            reverse_depth: 0,
            any_valid: false,
        }
    }
}

impl RawPosition {
    fn record_seen(&mut self, base: u8) {
        let index = match base {
            b'I' => 16,
            b'D' => 17,
            _ => base as usize,
        };
        if self.first_seen[index] == u32::MAX {
            self.first_seen[index] = self.next_order;
            self.next_order += 1;
        }
    }

    fn record_alt(&mut self, key: String) {
        if !self.alt_order.contains(&key) {
            self.alt_order.push(key);
        }
    }
}

fn reference_index(base: u8) -> usize {
    match base {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

fn finish_position(
    raw: RawPosition,
    position: i64,
    ref_base: u8,
    options: &PileupOptions,
    previous_position: &mut Option<i64>,
    contiguous: &mut usize,
) -> Option<PositionCounts> {
    if !raw.any_valid {
        return None;
    }
    let current_contiguous = match *previous_position {
        Some(previous) if previous + 1 == position => *contiguous + 1,
        _ => 0,
    };
    *previous_position = Some(position);
    *contiguous = current_contiguous;
    let RawPosition {
        mut matrix,
        deletions_f,
        deletions_r,
        insertions_f,
        insertions_r,
        insertions_all,
        alt_order,
        first_seen,
        depth,
        deleted_bases,
        forward_depth,
        reverse_depth,
        ..
    } = raw;
    let ref_index = reference_index(ref_base);
    for (_, count) in &deletions_f {
        matrix[6] += *count as i32;
    }
    if let Some(best) = deletions_f.iter().map(|(_, count)| count).max() {
        matrix[7] = *best as i32;
    }
    for (_, count) in &deletions_r {
        matrix[15] += *count as i32;
    }
    if let Some(best) = deletions_r.iter().map(|(_, count)| count).max() {
        matrix[16] = *best as i32;
    }
    let (ins_f_sum, ins_f_best) = insertion_stats(&insertions_f);
    let (ins_r_sum, ins_r_best) = insertion_stats(&insertions_r);
    matrix[4] = ins_f_sum as i32;
    matrix[5] = ins_f_best as i32;
    matrix[13] = ins_r_sum as i32;
    matrix[14] = ins_r_best as i32;

    let forward_sum: i64 = forward_depth + (0..4).map(|index| matrix[index] as i64).sum::<i64>();
    let reverse_sum: i64 =
        reverse_depth + (0..4).map(|index| matrix[index + 9] as i64).sum::<i64>();
    let mut alt_count = 0_i64;
    let mut major_alt_base = 0_u8;
    let mut all_alt_count = 0_u64;
    let del_count: i64 = deletions_f.iter().map(|(_, count)| *count).sum::<i64>()
        + deletions_r.iter().map(|(_, count)| *count).sum::<i64>();
    let ins_count = ins_f_sum + ins_r_sum;
    let mut pileup_list = [(0_u8, 0_i64); 6];
    let mut pileup_len = 0;
    for (index, base) in b"ACGT".iter().enumerate() {
        let current = matrix[index] as i64 + matrix[index + 9] as i64;
        if current > 0 {
            pileup_list[pileup_len] = (*base, current);
            pileup_len += 1;
        }
        if index != ref_index && current > alt_count {
            major_alt_base = *base;
            all_alt_count += current as u64;
            alt_count = current;
        }
    }
    let ref_count = depth
        - deleted_bases
        - (0..4)
            .filter(|index| *index != ref_index)
            .map(|index| matrix[index] as i64 + matrix[index + 9] as i64)
            .sum::<i64>();
    if ins_count > 0 {
        pileup_list[pileup_len] = (b'I', ins_count);
        pileup_len += 1;
    }
    if del_count > 0 {
        pileup_list[pileup_len] = (b'D', del_count);
        pileup_len += 1;
    }
    let pileup_list = &mut pileup_list[..pileup_len];
    pileup_list.sort_by(|left, right| {
        right.1.cmp(&left.1).then_with(|| {
            let order = |base: u8| {
                let index = match base {
                    b'I' => 16,
                    b'D' => 17,
                    _ => base as usize,
                };
                first_seen
                    .get(index)
                    .copied()
                    .filter(|order| *order != u32::MAX)
            };
            order(left.0).cmp(&order(right.0))
        })
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
        pass_af &= current_contiguous >= FLANKING_BASES;
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
    Some(PositionCounts {
        position,
        matrix,
        candidate,
        gvcf_ref: ref_count as u64,
        gvcf_total,
    })
}

fn deletion_bases(
    reference: &Reference,
    region: &Region,
    ref_seq: &[u8],
    position: i64,
    length: usize,
) -> Result<String> {
    let offset = (position - region.start) as usize;
    if offset + 1 + length <= ref_seq.len() {
        Ok(String::from_utf8_lossy(&ref_seq[offset + 1..offset + 1 + length]).into_owned())
    } else {
        Ok(String::from_utf8(reference.fetch_padded(
            &region.contig,
            position + 1,
            position + 1 + length as i64,
            b'A',
        )?)?)
    }
}

#[allow(clippy::too_many_arguments)]
fn add_indel(
    record: &Record,
    raw: &mut RawPosition,
    cigar: &[Cigar],
    index: usize,
    position: i64,
    query_position: usize,
    is_del: bool,
    reverse: bool,
    ref_base: u8,
    reference: &Reference,
    region: &Region,
    ref_seq: &[u8],
) -> Result<()> {
    let operation = cigar[index];
    let next = cigar.get(index + 1);
    if matches!(next, Some(Cigar::Del(_))) && !matches!(operation, Cigar::Del(_)) {
        let length: usize = cigar[index + 1..]
            .iter()
            .take_while(|operation| matches!(operation, Cigar::Del(_)))
            .map(|operation| operation.len() as usize)
            .sum();
        let target = if reverse {
            &mut raw.deletions_r
        } else {
            &mut raw.deletions_f
        };
        add_deletion(target, length);
        let deleted = deletion_bases(reference, region, ref_seq, position, length)?;
        raw.record_seen(b'D');
        raw.record_alt(format!("D{length}:{deleted}"));
    } else {
        let mut length = 0_usize;
        for operation in &cigar[index + 1..] {
            match operation {
                Cigar::Ins(size) => length += *size as usize,
                Cigar::Pad(_) => {}
                _ => break,
            }
        }
        if length > 0 {
            raw.record_seen(b'I');
            let first = usize::from(!is_del);
            let mut insertion = String::new();
            for offset in first..length + first {
                let source = query_position + offset;
                if source >= record.seq_len() {
                    bail!("insertion exceeds read length at {position}");
                }
                insertion.push(NT16[record.seq().encoded_base(source) as usize] as char);
            }
            let target = if reverse {
                &mut raw.insertions_r
            } else {
                &mut raw.insertions_f
            };
            add_insertion(target, &insertion);
            *raw.insertions_all.entry(insertion.clone()) += 1;
            raw.record_alt(format!("I{}{}", ref_base as char, insertion));
        }
    }
    Ok(())
}

fn add_matching_range(deltas: &mut [[i32; 2]], start: usize, end: usize, reverse: bool) {
    if start < end {
        let strand = usize::from(reverse);
        deltas[start][strand] += 1;
        deltas[end][strand] -= 1;
    }
}

pub fn create_direct_pileup(
    reader: &mut IndexedReader,
    reference: &Reference,
    region: &Region,
    options: &PileupOptions,
    pad_head: bool,
    pad_tail: bool,
) -> Result<PileupBatch> {
    let ref_seq = reference.fetch_region(region)?;
    reader.fetch((region.contig.as_bytes(), region.start, region.end))?;
    let mut pending: Vec<Record> = Vec::new();
    let mut next_record = None;
    let mut active = 0_usize;
    let mut active_ends: BTreeMap<i64, usize> = BTreeMap::new();
    let mut rows = Vec::new();
    let mut previous_position = None;
    let mut contiguous = 0;
    for window_start in (region.start..region.end).step_by(4096) {
        let window_end = (window_start + 4096).min(region.end);
        pending.retain(|record| record.cigar().end_pos() > window_start);
        loop {
            let record = if let Some(record) = next_record.take() {
                record
            } else {
                let mut record = Record::new();
                match reader.read(&mut record) {
                    Some(Ok(())) => record,
                    Some(Err(error)) => return Err(error.into()),
                    None => break,
                }
            };
            if record.pos() >= window_end {
                next_record = Some(record);
                break;
            }
            let flags = record.flags();
            if record.mapq() < options.min_mq
                || flags & (0x4 | 0x8 | 0x100 | 0x800) != 0
                || (flags & 0x1 != 0 && flags & 0x2 == 0)
            {
                continue;
            }
            while let Some((&end, _)) = active_ends.first_key_value()
                && end <= record.pos()
            {
                active -= active_ends.remove(&end).unwrap();
            }
            if active >= 1 << 20 {
                continue;
            }
            let end = record.cigar().end_pos();
            *active_ends.entry(end).or_default() += 1;
            active += 1;
            pending.push(record);
        }
        let mut raw = Vec::new();
        raw.resize_with((window_end - window_start) as usize, RawPosition::default);
        let mut matching_deltas = vec![[0_i32; 2]; (window_end - window_start) as usize + 1];
        for record in &pending {
            let operations = record.cigar().iter().copied().collect::<Vec<_>>();
            let reverse = record.flags() & 0x10 != 0;
            let sequence = record.seq();
            let mut ref_position = record.pos();
            let mut query_position = 0_usize;
            for (index, operation) in operations.iter().copied().enumerate() {
                let length = operation.len() as i64;
                match operation {
                    Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                        let overlap_start = ref_position.max(window_start);
                        let overlap_end = (ref_position + length).min(window_end);
                        let mut matching_start = overlap_start;
                        for position in overlap_start..overlap_end {
                            let query = query_position + (position - ref_position) as usize;
                            let encoded = sequence.encoded_base(query);
                            let ref_base = ref_seq[(position - region.start) as usize];
                            if !options.illumina
                                && NT16[encoded as usize] == ref_base
                                && matches!(ref_base, b'A' | b'C' | b'G' | b'T')
                            {
                                continue;
                            }
                            if !options.illumina && matching_start < position {
                                add_matching_range(
                                    &mut matching_deltas,
                                    (matching_start - window_start) as usize,
                                    (position - window_start) as usize,
                                    reverse,
                                );
                            }
                            matching_start = position + 1;
                            let row = &mut raw[(position - window_start) as usize];
                            row.any_valid = true;
                            if let Some(channel) = base_index(encoded, reverse) {
                                row.depth += 1;
                                row.matrix[channel] += 1;
                                let base = NT16[encoded as usize];
                                row.record_seen(encoded);
                                if base != ref_base {
                                    row.record_alt(format!("X{}", base as char));
                                }
                            }
                        }
                        if !options.illumina && matching_start < overlap_end {
                            add_matching_range(
                                &mut matching_deltas,
                                (matching_start - window_start) as usize,
                                (overlap_end - window_start) as usize,
                                reverse,
                            );
                        }
                        query_position += length as usize;
                    }
                    Cigar::Del(_) | Cigar::RefSkip(_) => {
                        let is_refskip = matches!(operation, Cigar::RefSkip(_));
                        for position in
                            ref_position.max(window_start)..(ref_position + length).min(window_end)
                        {
                            let row = &mut raw[(position - window_start) as usize];
                            row.any_valid = true;
                            if !is_refskip {
                                row.depth += 1;
                                row.deleted_bases += 1;
                                row.matrix[if reverse { 17 } else { 8 }] += 1;
                            }
                        }
                    }
                    Cigar::Ins(_) | Cigar::SoftClip(_) => {
                        query_position += length as usize;
                    }
                    Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
                if matches!(
                    operation,
                    Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_)
                ) {
                    continue;
                }
                let position = ref_position + length - 1;
                if position >= window_start
                    && position < window_end
                    && !matches!(operation, Cigar::RefSkip(_))
                {
                    let ref_base = ref_seq[(position - region.start) as usize].to_ascii_uppercase();
                    let query = if matches!(operation, Cigar::Del(_)) {
                        query_position
                    } else {
                        query_position - 1
                    };
                    add_indel(
                        record,
                        &mut raw[(position - window_start) as usize],
                        &operations,
                        index,
                        position,
                        query,
                        matches!(operation, Cigar::Del(_)),
                        reverse,
                        ref_base,
                        reference,
                        region,
                        &ref_seq,
                    )?;
                }
                ref_position += length;
            }
        }
        let mut matching_depth = [0_i64; 2];
        for (index, mut counts) in raw.into_iter().enumerate() {
            for strand in 0..2 {
                matching_depth[strand] += i64::from(matching_deltas[index][strand]);
            }
            let covered = matching_depth[0] + matching_depth[1];
            counts.forward_depth = matching_depth[0];
            counts.reverse_depth = matching_depth[1];
            counts.depth += covered;
            counts.any_valid |= covered > 0;
            let position = window_start + index as i64;
            let ref_base = ref_seq[(position - region.start) as usize].to_ascii_uppercase();
            if let Some(row) = finish_position(
                counts,
                position,
                ref_base,
                options,
                &mut previous_position,
                &mut contiguous,
            ) {
                rows.push(row);
            }
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
