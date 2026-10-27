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

use std::collections::{HashMap, VecDeque};

use super::ssw::smith_waterman;

const KMER_SIZE: usize = 32;
const MAX_MISMATCHES: i32 = 2;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RealignmentRead {
    pub sequence: Vec<u8>,
    pub position: i64,
    pub cigar: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Operation {
    Match,
    Insert,
    Delete,
    SoftClip,
    HardClip,
}

#[derive(Debug, Clone, Copy)]
struct CigarOp {
    operation: Operation,
    length: i32,
}

#[derive(Debug, Clone)]
struct ReadAlignment {
    position: i32,
    cigar: String,
    score: i32,
}

impl Default for ReadAlignment {
    fn default() -> Self {
        Self {
            position: -1,
            cigar: String::new(),
            score: 0,
        }
    }
}

#[derive(Debug)]
struct HaplotypeAlignment {
    score: i32,
    reads: Vec<ReadAlignment>,
    is_reference: bool,
    cigar: VecDeque<CigarOp>,
    reference_position: i64,
    position_map: Vec<i64>,
}

fn parse_cigar(cigar: &str) -> VecDeque<CigarOp> {
    let bytes = cigar.as_bytes();
    let mut result = VecDeque::new();
    let mut index = 0;
    while index < bytes.len() {
        if !bytes[index].is_ascii_digit() {
            index += 1;
            continue;
        }
        let start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == bytes.len() {
            break;
        }
        let operation = match bytes[index].to_ascii_uppercase() {
            b'=' | b'X' => Some(Operation::Match),
            b'I' => Some(Operation::Insert),
            b'D' => Some(Operation::Delete),
            b'S' => Some(Operation::SoftClip),
            b'H' => Some(Operation::HardClip),
            _ => None,
        };
        if let Some(operation) = operation {
            result.push_back(CigarOp {
                operation,
                length: cigar[start..index].parse().unwrap_or(0),
            });
        }
        index += 1;
    }
    result
}

fn cigar_string(cigar: &VecDeque<CigarOp>) -> String {
    cigar
        .iter()
        .map(|op| {
            let operation = match op.operation {
                Operation::Match => 'X',
                Operation::Insert => 'I',
                Operation::Delete => 'D',
                Operation::SoftClip => 'S',
                Operation::HardClip => 'H',
            };
            format!("{}{}", op.length, operation)
        })
        .collect()
}

fn fast_score(first: &[u8], second: &[u8]) -> (i32, i32) {
    let mut matches = 0;
    let mut mismatches = 0;
    for (left, right) in first.iter().zip(second) {
        if left != right && *left != b'N' && *right != b'N' {
            mismatches += 1;
            if mismatches == MAX_MISMATCHES + 1 {
                return (0, mismatches);
            }
        } else {
            matches += 1;
        }
    }
    (matches * 4 - mismatches * 6, mismatches)
}

fn fast_alignments(
    reads: &[RealignmentRead],
    reference: &[u8],
    haplotypes: &[Vec<u8>],
    prefix_length: usize,
    suffix_length: usize,
) -> Vec<HaplotypeAlignment> {
    let mut index: HashMap<Vec<u8>, Vec<(usize, usize)>> = HashMap::new();
    for (read_index, read) in reads.iter().enumerate() {
        if read.sequence.len() <= KMER_SIZE {
            continue;
        }
        for position in 0..=read.sequence.len() - KMER_SIZE {
            index
                .entry(read.sequence[position..position + KMER_SIZE].to_vec())
                .or_default()
                .push((read_index, position));
        }
    }
    let mut result = Vec::with_capacity(haplotypes.len());
    for haplotype in haplotypes {
        let is_reference = haplotype == reference;
        let mut coverage = vec![0_usize; haplotype.len()];
        let mut read_alignments = vec![ReadAlignment::default(); reads.len()];
        let mut haplotype_score = 0;
        if haplotype.len() >= KMER_SIZE {
            for position in 0..=haplotype.len() - KMER_SIZE {
                if let Some(occurrences) = index.get(&haplotype[position..position + KMER_SIZE]) {
                    for (read_index, read_position) in occurrences {
                        let target_start = position.saturating_sub(*read_position);
                        let read_length = reads[*read_index].sequence.len();
                        if target_start + read_length > haplotype.len() {
                            continue;
                        }
                        let alignment = &mut read_alignments[*read_index];
                        if alignment.position != -1 && alignment.position == target_start as i32 {
                            continue;
                        }
                        let (score, mismatches) = fast_score(
                            &haplotype[target_start..target_start + read_length],
                            &reads[*read_index].sequence,
                        );
                        if mismatches <= MAX_MISMATCHES {
                            let old_score = alignment.score;
                            coverage[target_start..target_start + read_length]
                                .iter_mut()
                                .for_each(|value| *value += 1);
                            if old_score < score {
                                alignment.score = score;
                                haplotype_score += score - old_score;
                                alignment.position = target_start as i32;
                                alignment.cigar = format!("{read_length}=");
                            }
                        }
                    }
                    if coverage[position] == 0
                        && position >= prefix_length
                        && position < haplotype.len().saturating_sub(suffix_length)
                        && !is_reference
                    {
                        haplotype_score = 0;
                        break;
                    }
                }
            }
        }
        if haplotype_score == 0 {
            read_alignments.fill(ReadAlignment::default());
        }
        result.push(HaplotypeAlignment {
            score: haplotype_score,
            reads: read_alignments,
            is_reference: false,
            cigar: VecDeque::new(),
            reference_position: 0,
            position_map: Vec::new(),
        });
    }
    result
}

fn position_map(haplotype_length: usize, cigar: &VecDeque<CigarOp>) -> Vec<i64> {
    let mut result = vec![0; haplotype_length];
    let mut shift = 0;
    let mut position = 0;
    for op in cigar {
        match op.operation {
            Operation::Match => {
                for _ in 0..op.length {
                    if position < result.len() {
                        result[position] = shift;
                    }
                    position += 1;
                }
            }
            Operation::SoftClip => {
                shift -= op.length as i64;
                for _ in 0..op.length {
                    if position < result.len() {
                        result[position] = shift;
                    }
                    position += 1;
                }
            }
            Operation::Delete => shift += op.length as i64,
            Operation::Insert => {
                for _ in 0..op.length {
                    if position < result.len() {
                        result[position] = shift;
                    }
                    shift -= 1;
                    position += 1;
                }
            }
            Operation::HardClip => {}
        }
    }
    result
}

fn aligned_length(cigar: &VecDeque<CigarOp>) -> i32 {
    cigar
        .iter()
        .filter(|op| op.operation != Operation::Delete)
        .map(|op| op.length)
        .sum()
}

fn merge_cigar(op: CigarOp, read_length: i32, cigar: &mut VecDeque<CigarOp>) {
    let before = aligned_length(cigar);
    let length = if op.operation == Operation::Delete {
        op.length
    } else {
        op.length.min(read_length - before)
    };
    if length <= 0 || before == read_length {
        return;
    }
    if let Some(last) = cigar.back_mut()
        && last.operation == op.operation
    {
        last.length += length;
        return;
    }
    cigar.push_back(CigarOp {
        operation: op.operation,
        length,
    });
}

fn push_front(op: CigarOp, cigar: &mut VecDeque<CigarOp>) {
    if op.length > 0 {
        cigar.push_front(op);
    }
}

fn left_trim(mut cigar: VecDeque<CigarOp>, read_position: i32) -> VecDeque<CigarOp> {
    let mut position = 0;
    while position != read_position {
        let Some(op) = cigar.pop_front() else {
            return cigar;
        };
        if matches!(
            op.operation,
            Operation::Match | Operation::SoftClip | Operation::Insert
        ) {
            if op.length + position > read_position {
                cigar.push_front(CigarOp {
                    operation: op.operation,
                    length: op.length - (read_position - position),
                });
            }
            position = (op.length + position).min(read_position);
        }
    }
    if cigar
        .front()
        .is_some_and(|op| op.operation == Operation::Delete)
    {
        cigar.pop_front();
    }
    cigar
}

fn is_match(op: CigarOp) -> bool {
    matches!(op.operation, Operation::Match | Operation::SoftClip)
}

fn composed_cigar(
    read_length: usize,
    read_alignment: &ReadAlignment,
    haplotype_cigar: &VecDeque<CigarOp>,
) -> VecDeque<CigarOp> {
    let read_length = read_length as i32;
    let mut read_cigar = parse_cigar(&read_alignment.cigar);
    let mut haplotype_cigar = left_trim(haplotype_cigar.clone(), read_alignment.position);
    let mut output = VecDeque::new();
    if read_cigar
        .front()
        .is_some_and(|op| op.operation == Operation::SoftClip)
    {
        merge_cigar(read_cigar[0], read_length, &mut output);
        read_cigar.pop_front();
    }
    while (!read_cigar.is_empty() || !haplotype_cigar.is_empty())
        && aligned_length(&output) < read_length
    {
        if !read_cigar.is_empty() && haplotype_cigar.is_empty() {
            merge_cigar(read_cigar[0], read_length, &mut output);
            read_cigar.pop_front();
            continue;
        }
        if read_cigar.is_empty() {
            break;
        }
        let mut read_op = read_cigar.pop_front().unwrap();
        let mut haplotype_op = haplotype_cigar.pop_front().unwrap();
        if is_match(read_op) && is_match(haplotype_op) {
            let length = read_op.length.min(haplotype_op.length);
            let operation = if read_op.operation == Operation::SoftClip
                || haplotype_op.operation == Operation::SoftClip
            {
                Operation::SoftClip
            } else {
                Operation::Match
            };
            merge_cigar(CigarOp { operation, length }, read_length, &mut output);
            read_op.length -= length;
            haplotype_op.length -= length;
            push_front(read_op, &mut read_cigar);
            push_front(haplotype_op, &mut haplotype_cigar);
        } else if read_op.operation == Operation::Delete && is_match(haplotype_op) {
            merge_cigar(read_op, read_length, &mut output);
            haplotype_op.length -= read_op.length;
            push_front(haplotype_op, &mut haplotype_cigar);
        } else if haplotype_op.operation == Operation::Delete && is_match(read_op) {
            merge_cigar(haplotype_op, read_length, &mut output);
            push_front(read_op, &mut read_cigar);
        } else if read_op.operation == Operation::Delete
            && haplotype_op.operation == Operation::Delete
        {
            merge_cigar(
                CigarOp {
                    operation: Operation::Delete,
                    length: read_op.length + haplotype_op.length,
                },
                read_length,
                &mut output,
            );
        } else if read_op.operation == Operation::Insert && is_match(haplotype_op) {
            read_op.length = (read_length - aligned_length(&output)).min(read_op.length);
            merge_cigar(read_op, read_length, &mut output);
            push_front(haplotype_op, &mut haplotype_cigar);
        } else if haplotype_op.operation == Operation::Insert && is_match(read_op) {
            haplotype_op.length = (read_length - aligned_length(&output)).min(haplotype_op.length);
            merge_cigar(haplotype_op, read_length, &mut output);
            read_op.length = (read_op.length - haplotype_op.length).max(0);
            push_front(read_op, &mut read_cigar);
        } else if read_op.operation == Operation::Insert
            && haplotype_op.operation == Operation::Insert
        {
            haplotype_op.length += read_op.length;
            merge_cigar(haplotype_op, read_length, &mut output);
        } else {
            return VecDeque::new();
        }
    }
    output
}

pub fn realign_reads(
    reads: &[RealignmentRead],
    reference: &[u8],
    haplotypes: &[Vec<u8>],
    reference_start: i64,
    prefix_length: usize,
    suffix_length: usize,
) -> Vec<RealignmentRead> {
    let mut alignments =
        fast_alignments(reads, reference, haplotypes, prefix_length, suffix_length);
    for (haplotype, alignment) in haplotypes.iter().zip(&mut alignments) {
        let haplotype_alignment = smith_waterman(haplotype, reference);
        if haplotype_alignment.score > 0 {
            alignment.is_reference = haplotype_alignment.cigar == format!("{}=", haplotype.len());
            alignment.cigar = parse_cigar(&haplotype_alignment.cigar);
            alignment.reference_position = haplotype_alignment.reference_start as i64;
        }
        alignment.position_map = position_map(haplotype.len(), &alignment.cigar);
    }
    for (read_index, read) in reads.iter().enumerate() {
        if alignments
            .iter()
            .any(|alignment| alignment.reads[read_index].score > 0)
        {
            continue;
        }
        for (haplotype, alignment) in haplotypes.iter().zip(&mut alignments) {
            if alignment.score == 0 {
                continue;
            }
            let read_alignment = smith_waterman(&read.sequence, haplotype);
            if read_alignment.score >= 1 && alignment.reads[read_index].score < read_alignment.score
            {
                alignment.reads[read_index] = ReadAlignment {
                    score: read_alignment.score,
                    position: read_alignment.reference_start,
                    cigar: read_alignment.cigar,
                };
            }
        }
    }
    alignments.sort_by_key(|alignment| alignment.score);
    let mut result = reads.to_vec();
    for (read_index, output) in result.iter_mut().enumerate() {
        let mut best_score = 0;
        let mut best_haplotype = None;
        for (haplotype_index, alignment) in alignments.iter().enumerate() {
            let score = alignment.reads[read_index].score;
            if score > best_score
                || (best_score > 0 && score == best_score && !alignment.is_reference)
            {
                best_score = score;
                best_haplotype = Some(haplotype_index);
            }
        }
        if let Some(haplotype_index) = best_haplotype {
            let alignment = &alignments[haplotype_index];
            let read_alignment = &alignment.reads[read_index];
            let position = read_alignment.position.max(0) as usize;
            if let Some(shift) = alignment.position_map.get(position) {
                let cigar = composed_cigar(
                    reads[read_index].sequence.len(),
                    read_alignment,
                    &alignment.cigar,
                );
                if !cigar.is_empty() {
                    output.position = reference_start
                        + alignment.reference_position
                        + read_alignment.position as i64
                        + shift;
                    output.cigar = cigar_string(&cigar);
                }
            }
        }
    }
    result
}
