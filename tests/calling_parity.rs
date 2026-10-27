use std::{
    fs,
    path::{Path, PathBuf},
};

use clair3::{
    call::{
        caller::{TensorInput, call_tensor_batch, write_vcf},
        decode::DecodeOptions,
    },
    io::{
        interval::Region,
        reference::Reference,
        variant::{VcfReader, VcfRecord},
    },
    model::{Model, ModelKind},
    tensor::{
        full_alignment::{FullAlignmentOptions, create_full_alignment},
        pileup::{PileupOptions, create_pileup},
    },
};
use tch::Device;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn records(path: &Path) -> Vec<String> {
    let mut reader = VcfReader::open(path).unwrap();
    let mut records = Vec::new();
    while let Some(record) = reader.next_record().unwrap() {
        records.push(record.to_line());
    }
    records
}

fn lines(records: &[VcfRecord]) -> Vec<String> {
    records.iter().map(VcfRecord::to_line).collect()
}

fn assert_records_match(actual: &[VcfRecord], expected_path: &Path) {
    let actual = lines(actual);
    let expected = records(expected_path);
    let expected_count = expected.len();
    let mut maximum_quality_difference = 0.0_f64;
    for (index, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
        let mut actual_record = VcfRecord::from_line(actual).unwrap();
        let mut expected_record = VcfRecord::from_line(expected).unwrap();
        let quality_difference = (actual_record.quality.parse::<f64>().unwrap()
            - expected_record.quality.parse::<f64>().unwrap())
        .abs();
        maximum_quality_difference = maximum_quality_difference.max(quality_difference);
        let gq_difference = (actual_record
            .sample_field(0, 1)
            .unwrap()
            .parse::<i32>()
            .unwrap()
            - expected_record
                .sample_field(0, 1)
                .unwrap()
                .parse::<i32>()
                .unwrap())
        .abs();
        actual_record.quality.clear();
        expected_record.quality.clear();
        for record in [&mut actual_record, &mut expected_record] {
            record.samples = record
                .samples
                .iter()
                .map(|sample| {
                    let mut fields: Vec<_> = sample.split(':').collect();
                    fields[1] = "";
                    fields.join(":")
                })
                .collect();
        }
        assert!(
            actual_record == expected_record && gq_difference <= 1,
            "record {index}: actual {actual}, expected {expected}; quality difference {quality_difference}, GQ difference {gq_difference}"
        );
    }
    assert_eq!(actual.len(), expected_count, "VCF record count differs");
    assert!(
        maximum_quality_difference <= 0.5,
        "maximum quality difference {maximum_quality_difference}"
    );
}

#[test]
#[ignore = "requires official demo inputs and legacy baselines"]
fn official_illumina_calls_match() {
    let root = root();
    let data = root.join("target/test-data/illumina");
    let reference = Reference::open(&data.join("GRCh38_chr20.fa")).unwrap();
    let region = Region::new("chr20", 0, reference.contig_len("chr20") as i64).unwrap();
    let bam = data.join("HG003_chr20_demo.bam");
    let baseline = root.join("target/legacy/baseline/official_ilmn");
    let confident_regions = clair3::io::interval::IntervalMap::from_bed(
        &baseline.join("tmp/split_beds/chr20"),
        Some("chr20"),
    )
    .unwrap();
    let pileup = create_pileup(
        &bam,
        &reference,
        &region,
        &PileupOptions {
            min_indel_af: 0.08,
            illumina: true,
            confident_regions: Some(&confident_regions),
            ..Default::default()
        },
    )
    .unwrap();
    let expected_info =
        fs::read_to_string(root.join("tests/fixtures/tensors/official_ilmn_pileup.info")).unwrap();
    let mut expected_positions = Vec::new();
    let mut expected_alt = Vec::new();
    let mut expected_indices = Vec::new();
    for (index, line) in expected_info.lines().enumerate() {
        let (position, alt) = line.split_once('\t').unwrap();
        let (prefix, reference) = position.rsplit_once(':').unwrap();
        let normalized = format!("{prefix}:{}", reference.as_bytes()[16] as char);
        if expected_positions.last() != Some(&normalized) {
            expected_indices.push(index as i64);
            expected_positions.push(normalized);
            expected_alt.push(alt.to_string());
        }
    }
    let expected_tensor =
        tch::Tensor::read_npy(root.join("tests/fixtures/tensors/official_ilmn_pileup.npy"))
            .unwrap()
            .index_select(0, &tch::Tensor::from_slice(&expected_indices));
    let expected_values = Vec::<i8>::try_from(&expected_tensor.flatten(0, -1)).unwrap();
    let actual_values = pileup.values_i8();
    let differences: Vec<_> = actual_values
        .iter()
        .zip(&expected_values)
        .enumerate()
        .filter(|(_, (left, right))| left != right)
        .take(20)
        .collect();
    assert_eq!(pileup.positions, expected_positions);
    assert_eq!(pileup.alt_info, expected_alt);
    assert_eq!(
        pileup.shape,
        expected_tensor
            .size()
            .into_iter()
            .map(|value| value as usize)
            .collect::<Vec<_>>()
    );
    assert!(
        differences.is_empty(),
        "first tensor differences: {differences:?}"
    );
    let model = Model::new(
        ModelKind::Pileup,
        Device::Cpu,
        Some(18),
        false,
        Some(&root.join("target/test-data/models/ilmn/pileup.pt")),
    )
    .unwrap();
    let actual = call_tensor_batch(
        &model,
        TensorInput {
            shape: &pileup.shape,
            values: &pileup.values,
            positions: &pileup.positions,
            alt_info: &pileup.alt_info,
        },
        &DecodeOptions {
            pileup: true,
            show_reference: true,
            ..Default::default()
        },
        200,
    )
    .unwrap();
    assert_records_match(&actual, &baseline.join("pileup.vcf.gz"));

    let realigned = root.join("target/verify/realigned_official_phased.bam");
    clair3::realign::workflow::realign_bam(
        &baseline.join("tmp/phase_output/phase_bam/chr20.bam"),
        &reference,
        "chr20",
        Some(100000),
        Some(300000),
        &realigned,
        &clair3::realign::workflow::RealignmentOptions {
            candidate_bed: Some(baseline.join("tmp/full_alignment_output/candidate_bed/chr20.0_1")),
            ..Default::default()
        },
    )
    .unwrap();
    let full_region = Region::new("chr20", 99967, 301033).unwrap();
    let candidate_pileup = create_pileup(
        &realigned,
        &reference,
        &full_region,
        &PileupOptions {
            min_indel_af: 0.15,
            exclude_deletions_from_af_denominator: true,
            allow_empty_flanking: true,
            call_head_tail: true,
            ..Default::default()
        },
    )
    .unwrap();
    let candidates: Vec<i64> = candidate_pileup
        .positions
        .iter()
        .map(|position| position.rsplit(':').nth(1).unwrap().parse::<i64>().unwrap() - 1)
        .collect();
    let full = create_full_alignment(
        &realigned,
        &reference,
        &full_region,
        &candidates,
        &[],
        &FullAlignmentOptions {
            matrix_depth: 55,
            illumina: true,
            ..Default::default()
        },
    )
    .unwrap();
    let model = Model::new(
        ModelKind::FullAlignment,
        Device::Cpu,
        Some(8),
        true,
        Some(&root.join("target/test-data/models/ilmn/full_alignment.pt")),
    )
    .unwrap();
    let options = DecodeOptions {
        add_indel_length: true,
        show_reference: true,
        ..Default::default()
    };
    let values: Vec<i32> = full.values.iter().map(|value| i32::from(*value)).collect();
    let actual = call_tensor_batch(
        &model,
        TensorInput {
            shape: &full.shape,
            values: &values,
            positions: &full.positions,
            alt_info: &full.alt_info,
        },
        &options,
        200,
    )
    .unwrap();
    assert_records_match(
        &actual,
        &baseline.join("tmp/full_alignment_output/full_alignment_chr20.0_1.vcf"),
    );
    let output = root.join("target/verify/official_illumina_full.vcf");
    fs::create_dir_all(output.parent().unwrap()).unwrap();
    write_vcf(&output, &reference, &actual, "HG003", None, false).unwrap();
    assert_records_match(&actual, &output);
}

#[test]
#[ignore = "requires official demo inputs and legacy baselines"]
fn production_full_call_includes_flanking_candidates() {
    let root = root();
    let baseline = root.join("target/legacy/baseline/official_ilmn");
    let bam = root.join("target/verify/realigned_official_phased.bam");
    let mut config = clair3::pipeline::call::CallConfig::new(
        bam,
        root.join("target/test-data/illumina/GRCh38_chr20.fa"),
        root.join("target/test-data/models/ilmn/full_alignment.pt"),
        "chr20".to_string(),
        clair3::config::Platform::Ilmn,
    );
    config.full_regions = Some(baseline.join("tmp/full_alignment_output/candidate_bed/chr20.0_1"));
    config.need_realignment = false;
    let actual = clair3::pipeline::call::call_var_bam(&config, false)
        .unwrap()
        .records;
    assert_records_match(
        &actual,
        &baseline.join("tmp/full_alignment_output/full_alignment_chr20.0_1.vcf"),
    );
}

#[test]
fn high_coverage_pileup_is_normalized_before_inference() {
    let checkpoint = root().join("target/test-data/models/ont/pileup.pt");
    if !checkpoint.exists() {
        return;
    }
    let model = Model::new(
        ModelKind::Pileup,
        Device::Cpu,
        Some(18),
        false,
        Some(&checkpoint),
    )
    .unwrap();
    let reference = Reference::open(&root().join("tests/fixtures/data/toy.fa")).unwrap();
    let region = Region::new("chrToy", 0, reference.contig_len("chrToy") as i64).unwrap();
    let batch = create_pileup(
        &root().join("tests/fixtures/data/toy.bam"),
        &reference,
        &region,
        &PileupOptions::default(),
    )
    .unwrap();
    let index = 7;
    let sample_size = 33 * 18;
    let values = &batch.values[index * sample_size..(index + 1) * sample_size];
    let input = (tch::Tensor::from_slice(values).reshape([1, 33, 18]) / (300.0 / 144.0))
        .to_kind(tch::Kind::Int);
    let expected = tch::Tensor::read_npy(
        root().join("tests/fixtures/models/npy/pileup_high_coverage_cpu.npy"),
    )
    .unwrap();
    let difference = (model.forward(&input, false) - expected)
        .abs()
        .max()
        .double_value(&[]);
    assert!(
        difference <= 1e-4,
        "model probabilities differ by {difference}"
    );
    let info =
        fs::read_to_string(root().join("tests/fixtures/tensors/pileup_default.info")).unwrap();
    let (position, alt) = info.lines().nth(7).unwrap().split_once('\t').unwrap();
    let positions = vec![position.to_string()];
    let alt_info = vec![alt.to_string()];
    let shape = vec![1, 33, 18];
    let actual = call_tensor_batch(
        &model,
        TensorInput {
            shape: &shape,
            values,
            positions: &positions,
            alt_info: &alt_info,
        },
        &DecodeOptions {
            pileup: true,
            ..Default::default()
        },
        1,
    )
    .unwrap();
    assert_eq!(actual.len(), 1);
    assert_eq!(actual[0].position, 2500);
    assert_eq!(actual[0].samples[0].split(':').next(), Some("0/1"));
    assert_eq!(actual[0].quality, "19.72");
}

#[test]
fn call_settings_reach_model_and_decoder() {
    let root = root();
    let checkpoint = root.join("target/test-data/models/ont/pileup.pt");
    let indel_checkpoint = root.join("tests/fixtures/models/pileup_indel.pt");
    if !checkpoint.exists() || !indel_checkpoint.exists() {
        return;
    }
    let data = root.join("tests/fixtures/data");
    let mut config = clair3::pipeline::call::CallConfig::new(
        data.join("toy.bam"),
        data.join("toy.fa"),
        checkpoint,
        "chrToy".to_string(),
        clair3::config::Platform::Ont,
    );
    config.start = Some(2484);
    config.end = Some(2516);
    config.quality_for_pass = Some(0.0);
    config.add_indel_length = false;
    config.haploid_sensitive = true;
    let sensitive = clair3::pipeline::call::call_var_bam(&config, true).unwrap();
    let variant = sensitive
        .records
        .iter()
        .find(|record| record.position == 2500)
        .unwrap();
    assert_eq!(variant.samples[0].split(':').next(), Some("1"));

    config.haploid_sensitive = false;
    config.haploid_precise = true;
    let precise = clair3::pipeline::call::call_var_bam(&config, true).unwrap();
    assert!(precise.records.iter().all(|record| record.position != 2500));

    config.haploid_precise = false;
    config.add_indel_length = true;
    config.checkpoint = indel_checkpoint;
    clair3::pipeline::call::call_var_bam(&config, true).unwrap();
}

#[test]
fn npy_metadata_matches_legacy_and_expands_for_calling() {
    let root = root();
    let checkpoint = root.join("target/test-data/models/ont/pileup.pt");
    if !checkpoint.exists() {
        return;
    }
    let data = root.join("tests/fixtures/data");
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let region = Region::new("chrToy", 0, reference.contig_len("chrToy") as i64).unwrap();
    let batch = create_pileup(
        &data.join("toy.bam"),
        &reference,
        &region,
        &PileupOptions::default(),
    )
    .unwrap();
    let tensor = clair3::pipeline::tensor::TensorData {
        shape: vec![1, 33, 18],
        values: batch.values[..33 * 18].to_vec(),
        positions: vec![batch.positions[0].clone()],
        alt_info: vec![batch.alt_info[0].clone()],
        nonvariants: Vec::new(),
    };
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("tensor.npy");
    tensor.write_npy(&path, &reference).unwrap();
    assert_eq!(
        fs::read_to_string(directory.path().join("tensor.npy.info")).unwrap(),
        format!("{}\t{}\n", batch.positions[0], batch.alt_info[0])
    );
    let mut restored = clair3::pipeline::tensor::TensorData::from_npy(&path).unwrap();
    restored.positions = restored.positions_with_reference(&reference).unwrap();
    assert_eq!(
        restored.positions,
        tensor.positions_with_reference(&reference).unwrap()
    );
    let records = clair3::pipeline::tensor::call_tensor_data(
        &checkpoint,
        ModelKind::Pileup,
        &restored,
        &DecodeOptions {
            pileup: true,
            show_reference: true,
            ..Default::default()
        },
        Device::Cpu,
    )
    .unwrap();
    assert!(!records.is_empty());
}

#[test]
fn native_min_bq_compatibility_is_stable() {
    let root = root();
    let data = root.join("tests/fixtures/data");
    let mut config = clair3::pipeline::call::CallConfig::new(
        data.join("toy.bam"),
        data.join("toy.fa"),
        PathBuf::new(),
        "chrToy".to_string(),
        clair3::config::Platform::Ont,
    );
    config.full_regions = Some(data.join("full_regions.bed"));
    config.no_phasing_for_fa = true;
    let full = clair3::pipeline::tensor::create_tensor_data(&config, false).unwrap();
    let pileup = clair3::pipeline::tensor::create_tensor_data(&config, true).unwrap();
    config.min_bq = 40;
    let filtered_full = clair3::pipeline::tensor::create_tensor_data(&config, false).unwrap();
    let filtered_pileup = clair3::pipeline::tensor::create_tensor_data(&config, true).unwrap();
    for (baseline, filtered) in [(&full, &filtered_full), (&pileup, &filtered_pileup)] {
        assert_eq!(baseline.values, filtered.values);
        assert_eq!(baseline.positions, filtered.positions);
        assert_eq!(baseline.alt_info, filtered.alt_info);
    }
}
