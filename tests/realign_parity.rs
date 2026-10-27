use std::{collections::BTreeSet, fs, path::PathBuf};

use clair3::realign::debruijn::consensus;
use serde::Deserialize;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[derive(Deserialize)]
struct ConsensusCase {
    reference: String,
    reads: String,
    base_quality: String,
    output: Vec<String>,
}

#[test]
fn consensus_matches_native_cases() {
    let text = fs::read_to_string(root().join("tests/fixtures/realign/consensus.json")).unwrap();
    let cases: Vec<ConsensusCase> = serde_json::from_str(&text).unwrap();
    assert_eq!(cases.len(), 238);
    for (case_index, case) in cases.iter().enumerate() {
        let reads: Vec<&[u8]> = case.reads.split(',').map(str::as_bytes).collect();
        let low_quality: Vec<BTreeSet<usize>> = case
            .base_quality
            .split(',')
            .map(|qualities| {
                qualities
                    .split_ascii_whitespace()
                    .map(|position| position.parse().unwrap())
                    .collect()
            })
            .collect();
        let actual = consensus(case.reference.as_bytes(), &reads, &low_quality);
        let expected: Vec<Vec<u8>> = case
            .output
            .iter()
            .map(|value| value.as_bytes().to_vec())
            .collect();
        if actual != expected {
            eprintln!(
                "case {case_index}: {} haplotypes, expected {}",
                actual.len(),
                expected.len()
            );
        }
        assert!(actual == expected, "consensus case {case_index}");
    }
}
