use std::{fs, path::PathBuf};

use clair3::{
    io::reference::Reference,
    realign::workflow::{RealignmentOptions, realign_bam},
};
use rust_htslib::bam::{self, Read as BamRead, Record};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[allow(clippy::type_complexity)]
fn snapshot(
    record: Record,
) -> (
    Vec<u8>,
    u16,
    i32,
    i64,
    u8,
    String,
    i32,
    i64,
    i64,
    Vec<u8>,
    Vec<u8>,
    Vec<(Vec<u8>, String)>,
) {
    let cigar = record
        .cigar()
        .iter()
        .map(|operation| operation.to_string())
        .collect();
    let auxiliary = record
        .aux_iter()
        .map(|item| {
            let (tag, value) = item.unwrap();
            (tag.to_vec(), format!("{value:?}"))
        })
        .collect();
    (
        record.qname().to_vec(),
        record.flags(),
        record.tid(),
        record.pos(),
        record.mapq(),
        cigar,
        record.mtid(),
        record.mpos(),
        record.insert_size(),
        record.seq().as_bytes(),
        record.qual().to_vec(),
        auxiliary,
    )
}

#[test]
#[ignore = "requires official demo data and native golden BAM"]
fn official_realigned_bam_matches() {
    let root = root();
    let input = root.join("target/test-data/illumina/HG003_chr20_demo.bam");
    let expected = root.join("tests/fixtures/realign/official_realigned.bam.100967_300033");
    assert!(input.exists() && expected.exists());
    let output = root.join("target/verify/realigned_official.bam");
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    let reference =
        Reference::open(&root.join("target/test-data/illumina/GRCh38_chr20.fa")).unwrap();
    let options = RealignmentOptions {
        candidate_bed: Some(root.join(
            "target/legacy/baseline/official_ilmn/tmp/full_alignment_output/candidate_bed/chr20.0_1",
        )),
        ..Default::default()
    };
    realign_bam(
        &input,
        &reference,
        "chr20",
        Some(100000),
        Some(300000),
        &output,
        &options,
    )
    .unwrap();
    assert!(output.with_extension("bam.bai").exists());
    let mut actual_reader = bam::Reader::from_path(&output).unwrap();
    let mut expected_reader = bam::Reader::from_path(&expected).unwrap();
    let actual: Vec<_> = actual_reader
        .records()
        .map(|record| snapshot(record.unwrap()))
        .collect();
    let expected: Vec<_> = expected_reader
        .records()
        .map(|record| snapshot(record.unwrap()))
        .collect();
    assert_eq!(actual.len(), expected.len());
    assert_eq!(actual.len(), 28730);
    for (index, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
        assert!(
            actual == expected,
            "BAM record {index}: actual {actual:?}; expected {expected:?}"
        );
    }
    let candidates = clair3::tensor::full_alignment::candidates_from_regions(
        options.candidate_bed.as_ref().unwrap(),
        "chr20",
    )
    .unwrap();
    let region =
        clair3::io::interval::Region::new("chr20", 0, reference.contig_len("chr20") as i64)
            .unwrap();
    let batch = clair3::tensor::full_alignment::create_full_alignment(
        &output,
        &reference,
        &region,
        &candidates,
        &[],
        &clair3::tensor::full_alignment::FullAlignmentOptions {
            matrix_depth: 55,
            illumina: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!candidates.is_empty());
    assert_eq!(batch.shape[0], candidates.len());
    assert_eq!(batch.positions.len(), candidates.len());
    assert!(batch.values.iter().any(|value| *value != 0));
}
