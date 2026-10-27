use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

use clair3::{
    io::{reference::Reference, variant::candidate_positions},
    tensor::pileup::{PileupBatch, PileupOptions, create_pileup},
};
use tch::Tensor;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn expected_info(case: &str) -> (Vec<String>, Vec<String>) {
    let text =
        fs::read_to_string(root().join(format!("tests/fixtures/tensors/{case}.info"))).unwrap();
    let mut positions = vec![];
    let mut alt_info = vec![];
    for line in text.lines() {
        let (position, alt) = line.split_once('\t').unwrap();
        positions.push(position.to_string());
        alt_info.push(alt.to_string());
    }
    (positions, alt_info)
}

fn expected_values(case: &str) -> (Vec<usize>, Vec<i8>) {
    let tensor =
        Tensor::read_npy(root().join(format!("tests/fixtures/tensors/{case}.npy"))).unwrap();
    let shape = tensor
        .size()
        .into_iter()
        .map(|value| value as usize)
        .collect();
    let values = Vec::<i8>::try_from(&tensor.flatten(0, -1)).unwrap();
    (shape, values)
}

fn check(batch: &PileupBatch, case: &str) {
    let (expected_shape, expected_values) = expected_values(case);
    assert_eq!(batch.shape, expected_shape, "shape for {case}");
    assert_eq!(
        batch.values_i8(),
        expected_values,
        "tensor values for {case}"
    );
    let (positions, alt_info) = expected_info(case);
    assert_eq!(batch.positions, positions, "positions for {case}");
    assert_eq!(batch.alt_info, alt_info, "alternate information for {case}");
}

fn generate(_case: &str, bam: &str, contig: &str, options: PileupOptions) -> PileupBatch {
    let reference = Reference::open(&root().join("tests/fixtures/data/toy.fa")).unwrap();
    let region =
        clair3::io::interval::Region::new(contig, 0, reference.contig_len(contig) as i64).unwrap();
    create_pileup(
        &root().join(format!("tests/fixtures/data/{bam}")),
        &reference,
        &region,
        &options,
    )
    .unwrap()
}

#[test]
fn default_bam_matches_golden() {
    check(
        &generate(
            "pileup_default",
            "toy.bam",
            "chrToy",
            PileupOptions::default(),
        ),
        "pileup_default",
    );
}

#[test]
fn head_tail_matches_golden() {
    let options = PileupOptions {
        call_head_tail: true,
        ..Default::default()
    };
    check(
        &generate("pileup_head_tail", "toy.bam", "chrToy", options),
        "pileup_head_tail",
    );
}

#[test]
fn long_indels_match_golden() {
    let options = PileupOptions {
        max_indel_length: 100_000,
        ..Default::default()
    };
    check(
        &generate("pileup_long", "toy.bam", "chrToy", options),
        "pileup_long",
    );
}

#[test]
fn known_sites_match_golden() {
    let path = root().join("tests/fixtures/data/known.vcf.gz");
    let known: BTreeSet<i64> = candidate_positions(&path, Some("chrToy"))
        .unwrap()
        .into_iter()
        .collect();
    let options = PileupOptions {
        max_indel_length: 100_000,
        known_positions: Some(&known),
        ..Default::default()
    };
    check(
        &generate("pileup_known", "toy.bam", "chrToy", options),
        "pileup_known",
    );
}

#[test]
fn cram_matches_golden() {
    check(
        &generate(
            "pileup_cram",
            "toy.cram",
            "chrToy",
            PileupOptions::default(),
        ),
        "pileup_cram",
    );
}

#[test]
fn colon_contig_matches_golden() {
    check(
        &generate(
            "pileup_colon",
            "toy.bam",
            "chr:colon",
            PileupOptions::default(),
        ),
        "pileup_colon",
    );
}

#[test]
fn hifi_matches_golden() {
    check(
        &generate("pileup_hifi", "toy.bam", "chrToy", PileupOptions::default()),
        "pileup_hifi",
    );
}

#[test]
fn empty_bam_has_no_candidates() {
    let batch = generate(
        "pileup_empty",
        "empty.bam",
        "chrEmpty",
        PileupOptions::default(),
    );
    assert!(batch.is_empty());
    assert_eq!(batch.shape[0], 0);
    assert!(!Path::new(&root().join("tests/fixtures/tensors/pileup_empty.npy")).exists());
}

#[test]
fn insertion_immediately_after_deletion_is_counted() {
    use rust_htslib::bam::{
        self, Header, Writer,
        record::{Cigar, CigarString, Record},
    };

    let directory = tempfile::tempdir().unwrap();
    let bam_path = directory.path().join("deletion_insertion.bam");
    let mut header = Header::new();
    let mut contig = bam::header::HeaderRecord::new(b"SQ");
    contig.push_tag(b"SN", "chrToy").push_tag(b"LN", 4096);
    header.push_record(&contig);
    let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    let mut record = Record::new();
    record.set(
        b"deletion_insertion",
        Some(&CigarString(vec![
            Cigar::Match(45),
            Cigar::Del(1),
            Cigar::Ins(3),
            Cigar::Match(45),
        ])),
        &[b'A'; 93],
        &[40; 93],
    );
    record.set_tid(0);
    record.set_pos(440);
    record.set_mapq(60);
    record.set_flags(0);
    record.set_mtid(-1);
    record.set_mpos(-1);
    record.set_insert_size(0);
    writer.write(&record).unwrap();
    drop(writer);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 0).unwrap();

    let data = root().join("tests/fixtures/data");
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let batch = create_pileup(
        &bam_path,
        &reference,
        &clair3::io::interval::Region::new("chrToy", 0, 4096).unwrap(),
        &PileupOptions {
            min_depth: 1,
            min_snp_af: 0.0,
            min_indel_af: 0.0,
            call_head_tail: true,
            ..Default::default()
        },
    )
    .unwrap();
    let anchor_base = reference.fetch("chrToy", 485, 486).unwrap()[0] as char;
    let index = batch
        .positions
        .iter()
        .position(|position| position == &format!("chrToy:486:{anchor_base}"))
        .unwrap();
    let expected = format!("I{anchor_base}AAA 1");
    assert!(batch.alt_info[index].contains(&expected));
}
