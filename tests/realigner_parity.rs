use std::{fs, path::PathBuf};

use clair3::realign::realigner::{RealignmentRead, realign_reads};
use serde::Deserialize;

#[derive(Deserialize)]
struct RealignmentCase {
    seqs: Vec<String>,
    positions: Vec<i64>,
    cigars: Vec<String>,
    reference: String,
    haplotypes: String,
    ref_start: i64,
    ref_prefix: usize,
    ref_suffix: usize,
    result_positions: Vec<i64>,
    result_cigars: Vec<String>,
}

#[test]
fn realignments_match_native_cases() {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/realign/realign_calls.json");
    let cases: Vec<RealignmentCase> =
        serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(cases.len(), 215);
    for (case_index, case) in cases.iter().enumerate() {
        let reads: Vec<RealignmentRead> = case
            .seqs
            .iter()
            .zip(&case.positions)
            .zip(&case.cigars)
            .map(|((sequence, position), cigar)| RealignmentRead {
                sequence: sequence.as_bytes().to_vec(),
                position: *position,
                cigar: cigar.clone(),
            })
            .collect();
        let haplotypes: Vec<Vec<u8>> = case
            .haplotypes
            .split_ascii_whitespace()
            .map(|haplotype| haplotype.as_bytes().to_vec())
            .collect();
        let actual = realign_reads(
            &reads,
            case.reference.as_bytes(),
            &haplotypes,
            case.ref_start,
            case.ref_prefix,
            case.ref_suffix,
        );
        let positions: Vec<i64> = actual.iter().map(|read| read.position).collect();
        let cigars: Vec<String> = actual
            .iter()
            .map(|read| read.cigar.replace('X', "M"))
            .collect();
        let expected: Vec<String> = case
            .result_cigars
            .iter()
            .map(|cigar| cigar.replace('X', "M"))
            .collect();
        if positions != case.result_positions || cigars != expected {
            let read_index = (0..reads.len())
                .find(|index| {
                    positions[*index] != case.result_positions[*index]
                        || cigars[*index] != expected[*index]
                })
                .unwrap();
            eprintln!(
                "case {case_index}, read {read_index}: position {} != {}, cigar {} != {}; input {} {}",
                positions[read_index],
                case.result_positions[read_index],
                cigars[read_index],
                expected[read_index],
                reads[read_index].position,
                reads[read_index].cigar
            );
            panic!("native realignment differs");
        }
    }
}
