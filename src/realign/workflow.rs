/*
Copyright 2020 Google LLC.
Copyright 2021 The University of Hong Kong, Department of Computer Science

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use rust_htslib::bam::{
    self, Format, Header, IndexedReader, Read as BamRead, Record, Writer,
    record::{Cigar, CigarString},
};

use crate::io::{compression::read_text_lines, interval::IntervalMap, reference::Reference};

use super::{debruijn::consensus, realigner::RealignmentRead};

const CHUNK_SIZE: i64 = 5000;
const MAX_WINDOW_SIZE: i64 = 1000;
const REGION_EXPANSION: i64 = 20;
const MIN_WINDOWS_DISTANCE: i64 = 80;
const MIN_PILEUP_MQ: u8 = 20;
const MIN_GRAPH_MQ: u8 = 14;
const MIN_GRAPH_BQ: u8 = 20;
const READ_NAME_CHARS: &[u8] =
    b"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$%&()*+./:;<=>?[]^`{|}~";

#[derive(Debug, Clone)]
pub struct RealignmentOptions {
    pub candidate_bed: Option<PathBuf>,
    pub extend_bed: Option<PathBuf>,
    pub min_coverage: f64,
    pub min_mq: u8,
    pub chunk_id: Option<usize>,
    pub chunk_num: Option<usize>,
}

impl Default for RealignmentOptions {
    fn default() -> Self {
        Self {
            candidate_bed: None,
            extend_bed: None,
            min_coverage: 2.0,
            min_mq: 5,
            chunk_id: None,
            chunk_num: None,
        }
    }
}

struct WorkRead {
    record: Record,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
    original_cigar: String,
    best_cigar: String,
    best_position: i64,
    best_score: Option<usize>,
    read_end: i64,
}

impl WorkRead {
    fn new(record: Record, sequence: Vec<u8>, qualities: Vec<u8>, cigar: String) -> Self {
        let read_end = record.pos()
            + sequence.len() as i64
            + bam_cigar(&cigar)
                .map(|cigar| {
                    cigar
                        .iter()
                        .map(|operation| match operation {
                            Cigar::Del(length) => *length as i64,
                            _ => 0,
                        })
                        .sum::<i64>()
                })
                .unwrap_or(0);
        Self {
            best_position: record.pos(),
            best_cigar: cigar.clone(),
            original_cigar: cigar,
            record,
            sequence,
            qualities,
            best_score: None,
            read_end,
        }
    }

    fn update(&mut self, position: i64, cigar: &str) {
        let cigar = cigar.replace('X', "M");
        if cigar == self.original_cigar && position == self.record.pos() {
            return;
        }
        if self.best_score.is_some() && cigar == self.best_cigar && position == self.best_position {
            return;
        }
        let score = bam_cigar(&cigar)
            .map(|cigar| {
                cigar
                    .iter()
                    .map(|operation| match operation {
                        Cigar::Ins(length) | Cigar::Del(length) => *length as usize,
                        _ => 0,
                    })
                    .sum()
            })
            .unwrap_or(0);
        if self.best_score.is_none_or(|best| score >= best) {
            self.best_position = position;
            self.best_cigar = cigar;
            self.best_score = Some(score);
        }
    }
}

fn bam_cigar(text: &str) -> Result<CigarString> {
    let bytes = text.as_bytes();
    let mut index = 0;
    let mut result = Vec::new();
    while index < bytes.len() {
        let start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if start == index || index == bytes.len() {
            bail!("invalid CIGAR '{text}'");
        }
        let length = text[start..index].parse()?;
        let operation = match bytes[index] {
            b'M' => Cigar::Match(length),
            b'I' => Cigar::Ins(length),
            b'D' => Cigar::Del(length),
            b'N' => Cigar::RefSkip(length),
            b'S' => Cigar::SoftClip(length),
            b'H' => Cigar::HardClip(length),
            b'P' => Cigar::Pad(length),
            b'=' => Cigar::Equal(length),
            b'X' => Cigar::Diff(length),
            _ => bail!("unsupported CIGAR operation in '{text}'"),
        };
        result.push(operation);
        index += 1;
    }
    Ok(CigarString(result))
}

fn short_read_name(mut index: usize) -> [u8; 5] {
    let mut result = [b'0'; 5];
    for position in (0..5).rev() {
        result[position] = READ_NAME_CHARS[index % READ_NAME_CHARS.len()];
        index /= READ_NAME_CHARS.len();
    }
    result
}

fn too_many_soft_clips(record: &Record) -> bool {
    let mut soft_clips = 0_u64;
    let mut total = 0_u64;
    for operation in record.cigar().iter() {
        total += operation.len() as u64;
        if matches!(operation, Cigar::SoftClip(_)) {
            soft_clips += operation.len() as u64;
        }
    }
    1.0 - (soft_clips as f64) / ((total + 1) as f64) < 0.55
}

fn reference_base(sequence: &[u8], reference_start: i64, position: i64) -> u8 {
    sequence
        .get((position - reference_start) as usize)
        .copied()
        .unwrap_or(b'N')
}

fn update_pileup(
    read: &WorkRead,
    pileup: &mut BTreeMap<i64, usize>,
    reference: &[u8],
    reference_start: i64,
    chunk_start: i64,
    chunk_end: i64,
) {
    if read.record.mapq() < MIN_PILEUP_MQ {
        return;
    }
    let mut reference_position = read.record.pos();
    let mut query_position = 0_usize;
    for operation in read.record.cigar().iter() {
        let length = operation.len() as usize;
        match operation {
            Cigar::Equal(_) => {
                reference_position += length as i64;
                query_position += length;
            }
            Cigar::Match(_) | Cigar::Diff(_) => {
                for _ in 0..length {
                    if read.qualities.get(query_position).copied().unwrap_or(0) >= MIN_GRAPH_BQ {
                        let reference_base =
                            reference_base(reference, reference_start, reference_position);
                        if matches!(reference_base, b'A' | b'C' | b'G' | b'T')
                            && read.sequence.get(query_position).copied().unwrap_or(b'N')
                                != reference_base
                        {
                            *pileup.entry(reference_position).or_default() += 1;
                        }
                    }
                    reference_position += 1;
                    query_position += 1;
                }
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) => {
                let previous_base =
                    reference_base(reference, reference_start, reference_position - 1);
                let outside = reference_position < chunk_start - REGION_EXPANSION
                    || reference_position > chunk_end + REGION_EXPANSION;
                let good_quality = read
                    .qualities
                    .get(query_position..query_position + length)
                    .is_some_and(|qualities| {
                        qualities.iter().all(|quality| *quality >= MIN_GRAPH_BQ)
                    });
                if !outside && matches!(previous_base, b'A' | b'C' | b'G' | b'T') && good_quality {
                    for position in
                        reference_position - length as i64..reference_position + length as i64
                    {
                        *pileup.entry(position).or_default() += 1;
                    }
                }
                query_position += length;
            }
            Cigar::Del(_) => {
                let previous_base =
                    reference_base(reference, reference_start, reference_position - 1);
                let outside = reference_position < chunk_start - REGION_EXPANSION
                    || reference_position > chunk_end + REGION_EXPANSION;
                if !outside && matches!(previous_base, b'A' | b'C' | b'G' | b'T') {
                    for position in reference_position..reference_position + length as i64 {
                        *pileup.entry(position).or_default() += 1;
                    }
                }
                reference_position += length as i64;
            }
            _ => {}
        }
    }
}

fn reference_slice(
    reference: &[u8],
    reference_start: i64,
    start: i64,
    end: i64,
) -> Option<Vec<u8>> {
    let (start, end) = if start <= end {
        (start, end)
    } else {
        (end, start)
    };
    let start = usize::try_from(start - reference_start).ok()?;
    let end = usize::try_from(end - reference_start).ok()?;
    reference.get(start..end).map(<[u8]>::to_vec)
}

fn write_read(read: &WorkRead, writer: &mut Writer) -> Result<()> {
    let mut record = read.record.clone();
    record.set_pos(read.best_position);
    record.set_cigar(Some(&bam_cigar(&read.best_cigar)?));
    writer.write(&record)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_chunk(
    reads: &mut Vec<WorkRead>,
    pileup: &mut BTreeMap<i64, usize>,
    reference: &[u8],
    reference_start: i64,
    chunk_start: i64,
    chunk_end: i64,
    options: &RealignmentOptions,
    writer: &mut Writer,
) -> Result<()> {
    let candidates: Vec<i64> = pileup
        .iter()
        .filter(|(position, count)| {
            (**count as f64) >= options.min_coverage
                && **position >= chunk_start - REGION_EXPANSION - 1
                && **position < chunk_end + REGION_EXPANSION
        })
        .map(|(position, _)| *position)
        .collect();
    if reads.is_empty() || candidates.is_empty() {
        return Ok(());
    }
    for split_index in 0..(chunk_end - chunk_start) / MAX_WINDOW_SIZE {
        let split_start = chunk_start + split_index * MAX_WINDOW_SIZE - REGION_EXPANSION - 1;
        let split_end = split_start + MAX_WINDOW_SIZE + REGION_EXPANSION * 2 + 1;
        let split_candidates: Vec<i64> = candidates
            .iter()
            .copied()
            .filter(|position| *position >= split_start && *position < split_end)
            .collect();
        let mut windows = Vec::new();
        for position in split_candidates {
            if let Some((_, end)) = windows.last_mut()
                && position <= *end + MIN_WINDOWS_DISTANCE * 2
            {
                *end = position;
                continue;
            }
            windows.push((position, position));
        }
        let windows: Vec<(i64, i64)> = windows
            .into_iter()
            .map(|(start, end)| (start - MIN_WINDOWS_DISTANCE, end + MIN_WINDOWS_DISTANCE))
            .collect();
        let Some(max_window_end) = windows.iter().map(|(_, end)| *end).max() else {
            continue;
        };
        let mut assignments: Vec<((i64, i64), Vec<usize>)> = windows
            .into_iter()
            .map(|window| (window, Vec::new()))
            .collect();
        for (read_index, read) in reads.iter().enumerate() {
            if read.record.pos() > max_window_end {
                continue;
            }
            let mut best = None;
            let mut best_overlap = 0;
            for (window_index, (window, _)) in assignments.iter().enumerate() {
                let overlap =
                    (read.read_end.min(window.1) - read.record.pos().max(window.0)).max(0);
                if overlap > best_overlap {
                    best_overlap = overlap;
                    best = Some(window_index);
                }
            }
            if let Some(window_index) = best {
                assignments[window_index].1.push(read_index);
            }
        }
        for ((start, end), assigned) in assignments {
            if end - start > MAX_WINDOW_SIZE {
                continue;
            }
            let Some(window_reference) = reference_slice(reference, reference_start, start, end)
            else {
                continue;
            };
            let graph_reads: Vec<&WorkRead> = assigned
                .iter()
                .map(|index| &reads[*index])
                .filter(|read| {
                    read.record.mapq() >= MIN_GRAPH_MQ
                        && read.record.pos() <= end
                        && read.read_end >= start
                })
                .collect();
            let sequences: Vec<&[u8]> = graph_reads
                .iter()
                .map(|read| read.sequence.as_slice())
                .collect();
            let low_quality: Vec<BTreeSet<usize>> = graph_reads
                .iter()
                .map(|read| {
                    read.qualities
                        .iter()
                        .enumerate()
                        .filter(|(_, quality)| **quality < 15)
                        .map(|(position, _)| position)
                        .collect()
                })
                .collect();
            let haplotypes = consensus(&window_reference, &sequences, &low_quality);
            if haplotypes.is_empty()
                || (haplotypes.len() == 1 && haplotypes[0] == window_reference)
                || assigned.is_empty()
            {
                continue;
            }
            let minimum_read_start = assigned
                .iter()
                .map(|index| reads[*index].record.pos())
                .min()
                .unwrap();
            let maximum_read_end = assigned
                .iter()
                .map(|index| reads[*index].read_end)
                .max()
                .unwrap();
            let temporary_start = (minimum_read_start.min(start) - REGION_EXPANSION).max(0);
            let temporary_end = maximum_read_end.max(end) + REGION_EXPANSION;
            let (Some(prefix), Some(center), Some(suffix)) = (
                reference_slice(reference, reference_start, temporary_start, start),
                reference_slice(reference, reference_start, start, end),
                reference_slice(reference, reference_start, end, temporary_end),
            ) else {
                continue;
            };
            let mut extended_reference = prefix.clone();
            extended_reference.extend_from_slice(&center);
            extended_reference.extend_from_slice(&suffix);
            let extended_haplotypes: Vec<Vec<u8>> = haplotypes
                .into_iter()
                .map(|haplotype| {
                    let mut extended = prefix.clone();
                    extended.extend_from_slice(&haplotype);
                    extended.extend_from_slice(&suffix);
                    extended
                })
                .collect();
            let input: Vec<RealignmentRead> = assigned
                .iter()
                .take(1000)
                .map(|index| {
                    let read = &reads[*index];
                    RealignmentRead {
                        sequence: read.sequence.clone(),
                        position: read.record.pos(),
                        cigar: read.original_cigar.clone(),
                    }
                })
                .collect();
            let realigned = super::realigner::realign_reads(
                &input,
                &extended_reference,
                &extended_haplotypes,
                temporary_start,
                prefix.len(),
                suffix.len(),
            );
            for (read_index, realigned) in assigned.iter().zip(realigned) {
                let read = &mut reads[*read_index];
                if realigned.cigar.is_empty()
                    || (realigned.cigar == read.original_cigar
                        && realigned.position == read.record.pos())
                {
                    continue;
                }
                read.update(realigned.position, &realigned.cigar);
            }
        }
    }
    let threshold = chunk_start - REGION_EXPANSION - MAX_WINDOW_SIZE;
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|index| reads[*index].best_position);
    let mut remove = vec![false; reads.len()];
    for index in order {
        if reads[index].best_position < threshold {
            write_read(&reads[index], writer)?;
            remove[index] = true;
        }
    }
    let mut index = 0;
    reads.retain(|_| {
        let keep = !remove[index];
        index += 1;
        keep
    });
    pileup.retain(|position, _| *position >= threshold);
    Ok(())
}

fn bed_range(path: &Path, contig: &str) -> Result<Option<(i64, i64)>> {
    let mut start = i64::MAX;
    let mut end = 0;
    for line in read_text_lines(path)? {
        let line = line?;
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns.len() < 3 || columns[0] != contig {
            continue;
        }
        start = start.min(columns[1].parse::<i64>()? + 1);
        end = end.max(columns[2].parse::<i64>()? + 1);
    }
    Ok((start != i64::MAX).then_some((start, end)))
}

pub fn realign_bam(
    input: &Path,
    reference: &Reference,
    contig: &str,
    contig_start: Option<i64>,
    contig_end: Option<i64>,
    output: &Path,
    options: &RealignmentOptions,
) -> Result<()> {
    if input == output {
        bail!("realigned BAM output must differ from the input BAM");
    }
    if contig_start.is_some() != contig_end.is_some() {
        bail!("contig start and end must be provided together");
    }
    if options.chunk_id.is_some() != options.chunk_num.is_some() {
        bail!("chunk ID and chunk count must be provided together");
    }
    let contig_length = reference.contig_len(contig) as i64;
    let range = if let Some(path) = &options.candidate_bed {
        bed_range(path, contig)?.with_context(|| {
            format!("contig {contig} is absent from BED file {}", path.display())
        })?
    } else if let (Some(chunk_id), Some(chunk_num)) = (options.chunk_id, options.chunk_num) {
        if chunk_id == 0 || chunk_num == 0 || chunk_id > chunk_num {
            bail!("chunk ID must be between 1 and the chunk count");
        }
        let chunk_size = if contig_length % chunk_num as i64 == 0 {
            contig_length / chunk_num as i64
        } else {
            contig_length / chunk_num as i64 + 1
        };
        let start = chunk_size * (chunk_id - 1) as i64;
        (start, start + chunk_size)
    } else if let (Some(start), Some(end)) = (contig_start, contig_end) {
        (start, end)
    } else {
        (1, contig_length)
    };
    let (range_start, range_end) = range;
    let has_range =
        options.candidate_bed.is_some() || options.chunk_id.is_some() || contig_start.is_some();
    let (reads_start, reads_end, reference_start, reference_end) = if has_range {
        (
            (range_start - MAX_WINDOW_SIZE).max(1) - 1,
            (range_end + MAX_WINDOW_SIZE).min(contig_length),
            (range_start - 5000).max(1) - 1,
            (range_end + 5000).min(contig_length),
        )
    } else {
        (0, contig_length, 0, contig_length)
    };
    let reference_sequence = reference.fetch(contig, reference_start, reference_end)?;
    let region_filter = if let Some(path) = &options.candidate_bed {
        Some(IntervalMap::from_bed(path, Some(contig))?)
    } else if let Some(path) = &options.extend_bed {
        Some(IntervalMap::from_bed(path, Some(contig))?)
    } else {
        None
    };
    let mut reader = IndexedReader::from_path(input)?;
    reader.set_reference(reference.path())?;
    let header = Header::from_template(reader.header());
    reader.fetch((contig.as_bytes(), reads_start, reads_end))?;
    let mut writer = Writer::from_path(output, &header, Format::Bam)?;
    let mut reads = Vec::new();
    let mut pileup = BTreeMap::new();
    let mut chunk: Option<(i64, i64)> = None;
    let mut next_read_name = 0_usize;
    for record in reader.records() {
        let mut record = record?;
        if record.mapq() < options.min_mq {
            continue;
        }
        let end = record.cigar().end_pos();
        if region_filter
            .as_ref()
            .is_some_and(|filter| !filter.overlaps(contig, record.pos(), end))
        {
            continue;
        }
        let name = short_read_name(next_read_name);
        next_read_name = (next_read_name + 1) % READ_NAME_CHARS.len().pow(5);
        if let Some((_, chunk_end)) = chunk
            && record.pos() >= chunk_end + REGION_EXPANSION
        {
            let (chunk_start, chunk_end) = chunk.unwrap();
            process_chunk(
                &mut reads,
                &mut pileup,
                &reference_sequence,
                reference_start,
                chunk_start,
                chunk_end,
                options,
                &mut writer,
            )?;
            chunk = Some((chunk_start + CHUNK_SIZE, chunk_end + CHUNK_SIZE));
        }
        if chunk.is_none() {
            chunk = Some((record.pos(), record.pos() + CHUNK_SIZE));
        }
        if record.cigar_len() == 0 || too_many_soft_clips(&record) {
            continue;
        }
        record.set_qname(&name);
        let mate_position = record.mpos();
        record.set_insert_size(mate_position + 1);
        let tags: Vec<Vec<u8>> = record
            .aux_iter()
            .filter_map(|item| item.ok().map(|(tag, _)| tag.to_vec()))
            .collect();
        for tag in tags {
            if tag != b"HP" {
                record.remove_aux(&tag)?;
            }
        }
        let sequence = record.seq().as_bytes();
        let qualities = record.qual().to_vec();
        let cigar = record
            .cigar()
            .iter()
            .map(|operation| operation.to_string())
            .collect();
        let read = WorkRead::new(record, sequence, qualities, cigar);
        update_pileup(
            &read,
            &mut pileup,
            &reference_sequence,
            reference_start,
            chunk.unwrap().0,
            chunk.unwrap().1,
        );
        reads.push(read);
    }
    if let Some((chunk_start, chunk_end)) = chunk {
        process_chunk(
            &mut reads,
            &mut pileup,
            &reference_sequence,
            reference_start,
            chunk_start,
            chunk_end,
            options,
            &mut writer,
        )?;
    }
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|index| reads[*index].best_position);
    for index in order {
        write_read(&reads[index], &mut writer)?;
    }
    drop(writer);
    bam::index::build(output, None::<&Path>, bam::index::Type::Bai, 1)?;
    Ok(())
}
