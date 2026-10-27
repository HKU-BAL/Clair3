use std::{
    fs,
    path::{Path, PathBuf},
};

use clair3::{
    io::{interval::IntervalMap, reference::Reference, variant::VcfRecord},
    pipeline::{
        merge::{MergeOptions, merge_illumina_records, merge_records, read_vcf},
        select::{SelectOptions, full_quality, phase_quality, select_candidates, select_het_snp},
    },
};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn lines(records: &[VcfRecord]) -> Vec<String> {
    records.iter().map(VcfRecord::to_line).collect()
}

fn expected(path: &Path) -> Vec<String> {
    lines(&read_vcf(path).unwrap().records)
}

#[test]
#[ignore = "requires official baseline outputs"]
fn candidate_and_merge_outputs_match() {
    let root = root();
    let baseline = root.join("target/legacy/baseline/official_ilmn");
    let pileup_path = baseline.join("pileup.vcf.gz");
    let reference =
        Reference::open(&root.join("target/test-data/illumina/GRCh38_chr20.fa")).unwrap();
    let output = root.join("target/verify/pipeline_ilmn/candidate_bed");
    let files = select_candidates(
        &pileup_path,
        None,
        &output,
        "chr20",
        "ilmn",
        &SelectOptions {
            qual_path: Some(baseline.join("tmp/full_alignment_output/candidate_bed/qual")),
            ..Default::default()
        },
        Some(&reference),
    )
    .unwrap();
    assert_eq!(files.len(), 1);
    for name in ["chr20.0_1", "full_aln_regions_chr20"] {
        assert_eq!(
            fs::read_to_string(output.join(name)).unwrap(),
            fs::read_to_string(
                baseline
                    .join("tmp/full_alignment_output/candidate_bed")
                    .join(name)
            )
            .unwrap(),
            "candidate file {name}"
        );
    }
    assert_eq!(
        fs::read_to_string(output.join("FULL_ALN_FILE_chr20")).unwrap(),
        format!("{}\n", output.join("chr20.0_1").display())
    );
    let pileup = read_vcf(&pileup_path).unwrap();
    let (variant_qual, reference_qual) = full_quality(&pileup.records, 0.3, Some(0.3), false);
    assert_eq!((variant_qual, reference_qual), (20.86, 21.98));
    let phase_qual = phase_quality(&pileup.records, Some(0.7), 0.3);
    let selected = select_het_snp(&pileup.records, "chr20", 0.7, Some(phase_qual));
    assert_eq!(
        lines(&selected),
        expected(&baseline.join("tmp/phase_output/phase_vcf/chr20.vcf"))
    );

    let full = read_vcf(&baseline.join("full_alignment.vcf.gz")).unwrap();
    let regions = IntervalMap::from_bed(
        &baseline.join("tmp/full_alignment_output/candidate_bed/full_aln_regions_chr20"),
        Some("chr20"),
    )
    .unwrap();
    let merged = merge_illumina_records(
        &pileup.records,
        &full.records,
        &regions,
        "chr20",
        &MergeOptions::default(),
    );
    assert_eq!(
        lines(&merged),
        expected(&baseline.join("tmp/merge_output/merge_chr20.vcf"))
    );

    let hifi = root.join("target/legacy/baseline/official_hifi");
    let merged = merge_records(
        &read_vcf(&hifi.join("pileup.vcf.gz")).unwrap().records,
        &read_vcf(&hifi.join("full_alignment.vcf.gz"))
            .unwrap()
            .records,
        &MergeOptions::default(),
    );
    assert_eq!(
        lines(&merged),
        expected(&hifi.join("tmp/merge_output/merge_chr20.vcf"))
    );
}

#[test]
fn small_contigs_use_one_chunk() {
    let root = root();
    let output = tempfile::tempdir().unwrap();
    let mut options = clair3::pipeline::run::RunOptions::new(
        root.join("tests/fixtures/data/toy.bam"),
        root.join("tests/fixtures/data/toy.fa"),
        root.join("target/test-data/models/ont"),
        output.path().to_path_buf(),
        clair3::config::Platform::Ont,
    );
    options.threads = 2;
    options.bed = Some(root.join("tests/fixtures/data/regions.bed"));
    options.include_all_contigs = true;
    options.no_phasing_for_fa = true;
    options.gvcf = true;
    options.output_all_contigs_in_gvcf_header = true;
    clair3::pipeline::run::run_pipeline(&options).unwrap();
    assert_eq!(
        fs::read_to_string(output.path().join("tmp/CHUNK_LIST")).unwrap(),
        "chr:colon 1 1\nchrToy 1 1\n"
    );
    let gvcf = read_vcf(&output.path().join("merge_output.gvcf.gz")).unwrap();
    assert!(gvcf.records.iter().any(|record| record.contig == "chrToy"
        && record.position == 2061
        && record.info == "END=2299"));
    assert!(
        gvcf.records
            .iter()
            .any(|record| record.contig == "chr:colon"
                && record.position == 201
                && record.info == "END=1024")
    );
}

#[test]
fn dwell_tag_validation_reads_sequential_bam() {
    let root = root();
    let reference = root.join("tests/fixtures/data/toy.fa");
    clair3::pipeline::run::validate_dwell_tags(
        &root.join("tests/fixtures/data/toy_dwell.bam"),
        &reference,
    )
    .unwrap();
    assert!(
        clair3::pipeline::run::validate_dwell_tags(
            &root.join("tests/fixtures/data/toy.bam"),
            &reference
        )
        .is_err()
    );
}

#[cfg(unix)]
#[test]
fn environment_checks_tool_versions_without_python() {
    use std::os::unix::fs::PermissionsExt;

    let root = root();
    let directory = tempfile::tempdir().unwrap();
    let version_script = |name: &str, version: &str| {
        let path = directory.path().join(name);
        fs::write(&path, format!("#!/bin/sh\nprintf '%s\\n' '{version}'\n")).unwrap();
        let mut permissions = fs::metadata(&path).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&path, permissions).unwrap();
        path.to_string_lossy().into_owned()
    };
    let samtools = version_script("samtools", "samtools 1.15.1");
    let whatshap = version_script("whatshap", "1.10");
    clair3::pipeline::run::validate_tool_version(&samtools, "samtools", &[1, 15, 1]).unwrap();
    clair3::pipeline::run::validate_tool_version(&whatshap, "whatshap", &[1, 7]).unwrap();
    let old_samtools = version_script("old-samtools", "samtools 1.15");
    let mut options = clair3::pipeline::run::RunOptions::new(
        root.join("tests/fixtures/data/toy.bam"),
        root.join("tests/fixtures/data/toy.fa"),
        root.join("target/test-data/models/ont"),
        directory.path().join("output"),
        clair3::config::Platform::Ont,
    );
    options.samtools = old_samtools;
    options.no_phasing_for_fa = true;
    let error = clair3::pipeline::run::check_environment(&options).unwrap_err();
    assert!(error.to_string().contains("too old"), "{error}");
}

#[test]
#[ignore = "requires external phasing tools"]
fn phasing_supports_colon_contigs() {
    let root = root();
    let output = tempfile::tempdir().unwrap();
    let mut options = clair3::pipeline::run::RunOptions::new(
        root.join("tests/fixtures/data/toy.bam"),
        root.join("tests/fixtures/data/toy.fa"),
        root.join("target/test-data/models/ilmn"),
        output.path().to_path_buf(),
        clair3::config::Platform::Ilmn,
    );
    options.threads = 2;
    options.bed = Some(root.join("tests/fixtures/data/regions.bed"));
    options.include_all_contigs = true;
    options.use_whatshap_for_final_phasing = true;
    options.use_whatshap_for_final_haplotagging = true;
    for (name, path) in [
        ("samtools", "samtools"),
        ("whatshap", "whatshap"),
        ("bgzip", "bgzip"),
        ("tabix", "tabix"),
    ] {
        let path = root.join("target/legacy/env/bin").join(path);
        let path = path.to_string_lossy().into_owned();
        match name {
            "samtools" => options.samtools = path,
            "whatshap" => options.whatshap = path,
            "bgzip" => options.bgzip = path,
            _ => options.tabix = path,
        }
    }
    clair3::pipeline::run::run_pipeline(&options).unwrap();
    assert!(
        output
            .path()
            .join("tmp/phase_output/phase_bam/chr:colon.bam.bai")
            .exists()
    );
    let merged = read_vcf(&output.path().join("merge_output.vcf.gz")).unwrap();
    assert!(
        merged
            .records
            .iter()
            .any(|record| record.contig == "chr:colon")
    );
}

#[test]
fn bp_resolution_preserves_legacy_counts_and_end() {
    let root = root();
    let mut config = clair3::pipeline::call::CallConfig::new(
        root.join("tests/fixtures/data/toy.bam"),
        root.join("tests/fixtures/data/toy.fa"),
        PathBuf::new(),
        "chrToy".to_string(),
        clair3::config::Platform::Ont,
    );
    config.start = Some(1);
    config.end = Some(4096);
    config.gvcf = true;
    config.gvcf_options = clair3::pipeline::gvcf::GvcfOptions {
        bp_resolution: true,
        ..Default::default()
    };
    let data = clair3::pipeline::tensor::create_tensor_data(&config, true).unwrap();
    assert_eq!(data.nonvariants.len(), 4095);
    assert_eq!(data.nonvariants.first().unwrap().position, 1);
    assert_eq!(data.nonvariants.last().unwrap().position, 4095);
    let site = data
        .nonvariants
        .iter()
        .find(|record| record.position == 1100)
        .unwrap();
    assert_eq!(site.samples[0], "./.:0:8:95,0,95");
}
