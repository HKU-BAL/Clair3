use std::{fs, path::PathBuf};

use clair3::realign::ssw::smith_waterman;
use serde::Deserialize;

#[derive(Deserialize)]
struct AlignmentCase {
    query: String,
    reference: String,
    score: i32,
    ref_begin: i32,
    ref_end: i32,
    query_begin: i32,
    query_end: i32,
    cigar: String,
}

#[test]
fn alignments_match_native_cases() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/realign/ssw.json");
    let cases: Vec<AlignmentCase> =
        serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap();
    assert!(!cases.is_empty());
    let mut failures = 0;
    for (index, case) in cases.iter().enumerate() {
        let actual = smith_waterman(case.query.as_bytes(), case.reference.as_bytes());
        if actual.score != case.score
            || actual.reference_start != case.ref_begin
            || actual.reference_end != case.ref_end
            || actual.query_start != case.query_begin
            || actual.query_end != case.query_end
            || actual.cigar != case.cigar
        {
            eprintln!(
                "case {index}: actual {actual:?}; expected score={}, ref={}..{}, query={}..{}, cigar={}",
                case.score,
                case.ref_begin,
                case.ref_end,
                case.query_begin,
                case.query_end,
                case.cigar
            );
            failures += 1;
            if failures == 5 {
                break;
            }
        }
    }
    assert_eq!(
        failures,
        0,
        "{} of {} native alignments differ",
        failures,
        cases.len()
    );
}
