use std::{fs, path::PathBuf};

use clair3::{
    io::reference::Reference,
    tensor::full_alignment::{
        FullAlignmentBatch, FullAlignmentOptions, candidates_from_regions_cffi,
        create_full_alignment, phased_variants_from_vcf,
    },
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

fn check(batch: &FullAlignmentBatch, case: &str) {
    let expected =
        Tensor::read_npy(root().join(format!("tests/fixtures/tensors/{case}.npy"))).unwrap();
    let shape = expected
        .size()
        .into_iter()
        .map(|value| value as usize)
        .collect::<Vec<_>>();
    assert_eq!(batch.shape, shape, "shape for {case}");
    let expected_values = Vec::<i8>::try_from(&expected.flatten(0, -1)).unwrap();
    let differences: Vec<(usize, i8, i8)> = batch
        .values
        .iter()
        .zip(&expected_values)
        .enumerate()
        .filter(|(_, (actual, expected))| actual != expected)
        .take(12)
        .map(|(index, (actual, expected))| (index, *actual, *expected))
        .collect();
    let mismatch_count = batch
        .values
        .iter()
        .zip(&expected_values)
        .filter(|(actual, expected)| actual != expected)
        .count();
    assert!(
        differences.is_empty(),
        "tensor values for {case}: {mismatch_count} mismatches, first (index, actual, expected): {differences:?}"
    );
    let (positions, alt_info) = expected_info(case);
    assert_eq!(batch.positions, positions, "positions for {case}");
    assert_eq!(batch.alt_info, alt_info, "alternate information for {case}");
}

fn generate(
    bam: &str,
    contig: &str,
    options: FullAlignmentOptions,
    phased: bool,
) -> FullAlignmentBatch {
    let data = root().join("tests/fixtures/data");
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let region =
        clair3::io::interval::Region::new(contig, 0, reference.contig_len(contig) as i64).unwrap();
    let candidates = candidates_from_regions_cffi(&data.join("full_regions.bed"), contig).unwrap();
    let variants = if phased {
        phased_variants_from_vcf(&data.join("phased.vcf.gz"), contig).unwrap()
    } else {
        vec![]
    };
    create_full_alignment(
        &data.join(bam),
        &reference,
        &region,
        &candidates,
        &variants,
        &options,
    )
    .unwrap()
}

#[test]
fn full_alignment_sequence_matches_golden() {
    check(
        &generate("toy.bam", "chrToy", FullAlignmentOptions::default(), false),
        "full_default",
    );
    let options = FullAlignmentOptions {
        max_indel_length: 100_000,
        ..Default::default()
    };
    check(&generate("toy.bam", "chrToy", options, false), "full_long");
    let options = FullAlignmentOptions {
        need_haplotagging: true,
        ..Default::default()
    };
    check(&generate("toy.bam", "chrToy", options, true), "full_phased");
    let options = FullAlignmentOptions {
        enable_dwell_time: true,
        ..Default::default()
    };
    check(
        &generate("toy_dwell.bam", "chrToy", options, false),
        "full_dwell",
    );
    check(
        &generate("toy.cram", "chrToy", FullAlignmentOptions::default(), false),
        "full_cram",
    );
    check(
        &generate(
            "toy.bam",
            "chr:colon",
            FullAlignmentOptions::default(),
            false,
        ),
        "full_colon",
    );
    let options = FullAlignmentOptions {
        matrix_depth: 55,
        ..Default::default()
    };
    check(&generate("toy.bam", "chrToy", options, false), "full_hifi");
    let data = root().join("tests/fixtures/data");
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let batch = create_full_alignment(
        &data.join("toy.bam"),
        &reference,
        &clair3::io::interval::Region::new("chrToy", 4, 22).unwrap(),
        &[4],
        &[],
        &FullAlignmentOptions::default(),
    )
    .unwrap();
    assert_eq!(batch.positions, vec!["chrToy:5:A"]);
}

#[test]
fn sampling_is_repeatable_across_concurrent_calls() {
    let calls: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(|| {
                generate("toy.bam", "chrToy", FullAlignmentOptions::default(), false)
            })
        })
        .collect();
    let batches: Vec<_> = calls.into_iter().map(|call| call.join().unwrap()).collect();
    for batch in &batches[1..] {
        assert_eq!(batch.positions, batches[0].positions);
        assert_eq!(batch.alt_info, batches[0].alt_info);
        let differences = batch
            .values
            .iter()
            .zip(&batches[0].values)
            .filter(|(left, right)| left != right)
            .count();
        assert_eq!(differences, 0, "sampled tensors differ between calls");
    }
}

#[test]
#[ignore = "requires official HiFi inputs and CFFI baseline tensors"]
fn official_hifi_phased_context_matches() {
    let root = root();
    let baseline = root.join("target/legacy/baseline/official_hifi");
    let candidate_path = baseline.join("tmp/full_alignment_output/candidate_bed/chr20.0_1");
    let content = fs::read_to_string(&candidate_path).unwrap();
    let mut start = i64::MAX;
    let mut end = 0;
    for line in content.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        start = start.min(fields[1].parse::<i64>().unwrap());
        end = end.max(fields[2].parse::<i64>().unwrap());
    }
    let data = root.join("target/test-data/pacbio_hifi");
    let reference = Reference::open(&data.join("GRCh38_no_alt_chr20.fa")).unwrap();
    let batch = create_full_alignment(
        &data.join("HG003_chr20_demo.bam"),
        &reference,
        &clair3::io::interval::Region::new("chr20", start, end).unwrap(),
        &candidates_from_regions_cffi(&candidate_path, "chr20").unwrap(),
        &phased_variants_from_vcf(
            &baseline.join("tmp/phase_output/phase_vcf/phased_chr20.vcf.gz"),
            "chr20",
        )
        .unwrap(),
        &FullAlignmentOptions {
            matrix_depth: 55,
            need_haplotagging: true,
            ..Default::default()
        },
    )
    .unwrap();
    let expected = Tensor::read_npy(root.join("target/verify/hifi-cffi.npy")).unwrap();
    assert_eq!(
        batch.shape,
        expected
            .size()
            .iter()
            .map(|value| *value as usize)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        batch.values,
        Vec::<i8>::try_from(&expected.flatten(0, -1)).unwrap()
    );
    let expected_info = fs::read_to_string(root.join("target/verify/hifi-cffi.npy.info")).unwrap();
    assert_eq!(
        batch.alt_info,
        expected_info
            .lines()
            .map(|line| line.split('\t').nth(1).unwrap().to_string())
            .collect::<Vec<_>>()
    );
}
