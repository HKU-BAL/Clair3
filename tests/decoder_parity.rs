use std::{fs, path::PathBuf};

use clair3::call::decode::{DecodeOptions, Prediction, decode_variant};
use serde::Deserialize;

#[derive(Deserialize)]
struct DecoderCase {
    contig: String,
    position: i64,
    reference: String,
    alt_info: String,
    gt21: Vec<f32>,
    genotype: Vec<f32>,
    length1: Vec<f32>,
    length2: Vec<f32>,
    show_reference: bool,
    haploid_precise: bool,
    haploid_sensitive: bool,
    qual: Option<f64>,
    add_indel_length: bool,
    gvcf: bool,
    pileup: bool,
    enable_long_indel: bool,
    maximum_variant_length: usize,
    keep_iupac: bool,
    output: Option<String>,
}

#[test]
fn decoder_matches_legacy_cases() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/call/decoder.json");
    let cases: Vec<DecoderCase> = serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(cases.len(), 630);
    let mut failures = 0;
    for (index, case) in cases.iter().enumerate() {
        let prediction = Prediction {
            gt21: case.gt21.clone(),
            genotype: case.genotype.clone(),
            length1: case.length1.clone(),
            length2: case.length2.clone(),
        };
        let options = DecodeOptions {
            show_reference: case.show_reference,
            haploid_precise: case.haploid_precise,
            haploid_sensitive: case.haploid_sensitive,
            quality_for_pass: case.qual,
            add_indel_length: case.add_indel_length,
            gvcf: case.gvcf,
            pileup: case.pileup,
            enable_long_indel: case.enable_long_indel,
            maximum_variant_length: case.maximum_variant_length,
            keep_iupac: case.keep_iupac,
        };
        let actual = decode_variant(
            &case.contig,
            case.position,
            &case.reference,
            16,
            &case.alt_info,
            &prediction,
            &options,
        )
        .unwrap()
        .map(|record| record.to_line() + "\n");
        if actual != case.output {
            eprintln!(
                "case {index}: actual {actual:?}, expected {:?}",
                case.output
            );
            failures += 1;
            if failures == 5 {
                break;
            }
        }
    }
    assert_eq!(failures, 0);
}
