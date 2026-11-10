use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::PathBuf,
    process::Command,
};

use clair3::{
    config::Platform,
    io::reference::Reference,
    pipeline::{call::CallConfig, unification::candidate_details_from_bam},
};

fn env_path(name: &str, default: &str) -> PathBuf {
    std::env::var_os(name)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default))
}

/// One candidate row: the depth and the reads supporting every allele.
type Candidate = (u64, BTreeMap<String, BTreeSet<String>>);

fn parse_candidates(text: &str) -> BTreeMap<u64, Candidate> {
    let mut rows = BTreeMap::new();
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 3, "invalid candidate row: {line}");
        let position: u64 = fields[0]
            .split_whitespace()
            .nth(1)
            .unwrap()
            .parse()
            .unwrap_or_else(|error| panic!("invalid candidate position in {line}: {error}"));
        let depth: u64 = fields[1]
            .parse()
            .unwrap_or_else(|error| panic!("invalid candidate depth in {line}: {error}"));
        let alleles: BTreeMap<String, String> = serde_json::from_str(fields[2])
            .unwrap_or_else(|error| panic!("invalid candidate JSON in {line}: {error}"));
        let alleles = alleles
            .into_iter()
            .map(|(allele, reads)| {
                (
                    allele,
                    reads.split_whitespace().map(str::to_string).collect(),
                )
            })
            .collect();
        rows.insert(position, (depth, alleles));
    }
    rows
}

#[test]
fn native_candidate_discovery_matches_the_reference() {
    let python = env_path(
        "CLAIR3_ORIGINAL_PYTHON",
        "/Users/Aqua/Temp/Clair3/target/tmp/macos-validation-venv/bin/python",
    );
    let original = env_path(
        "CLAIR3_ORIGINAL_ROOT",
        "/Users/Aqua/Temp/Clair3/target/tmp/original-main",
    )
    .join("clair3.py");
    let bam = env_path(
        "CLAIR3_UNIFICATION_BAM",
        "/Users/Aqua/Temp/Clair3/target/tmp/training-original-ru-inputs/chr21.ht.bam",
    );
    let reference_path = env_path(
        "CLAIR3_UNIFICATION_REFERENCE",
        "/Users/Aqua/Temp/clairBenchmark/GRCh38_no_alt_analysis_set.fasta",
    );
    let bed = env_path(
        "CLAIR3_UNIFICATION_BED",
        "/Users/Aqua/Temp/Clair3/target/tmp/training-chr21-inputs/truth/region.bed",
    );
    let contig = std::env::var("CLAIR3_UNIFICATION_CONTIG").unwrap_or_else(|_| "chr21".to_string());
    let platform =
        std::env::var("CLAIR3_UNIFICATION_PLATFORM").unwrap_or_else(|_| "hifi".to_string());
    let required = [
        python.clone(),
        original.clone(),
        bam.clone(),
        reference_path.clone(),
        bed.clone(),
    ];
    let missing: Vec<_> = required.iter().filter(|path| !path.exists()).collect();
    if !missing.is_empty() {
        eprintln!(
            "SKIP native_candidate_discovery_matches_the_reference: missing {}",
            missing
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return;
    }

    let directory = tempfile::tempdir().unwrap();
    let extended = directory.path().join("extended.bed");
    let output = Command::new(&python)
        .arg(&original)
        .arg("SplitExtendBed")
        .arg("--bed_fn")
        .arg(&bed)
        .arg("--output_fn")
        .arg(&extended)
        .arg("--ctgName")
        .arg(&contig)
        .arg("--expand_region_size")
        .arg("33")
        .output()
        .expect("failed to run the reference SplitExtendBed");
    assert!(
        output.status.success(),
        "the reference SplitExtendBed failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let (start, end) = {
        let mut start = i64::MAX;
        let mut end = i64::MIN;
        for line in fs::read_to_string(&extended).unwrap().lines() {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() >= 3 {
                start = start.min(fields[1].parse::<i64>().unwrap());
                end = end.max(fields[2].parse::<i64>().unwrap());
            }
        }
        (start, end)
    };

    // Reference candidate rows with the native allele-frequency threshold.
    let reference_rows = directory.path().join("reference-candidates.txt");
    let status = Command::new(&python)
        .arg(&original)
        .arg("CreateTensorFullAlignment")
        .arg("--bam_fn")
        .arg(&bam)
        .arg("--ref_fn")
        .arg(&reference_path)
        .arg("--ctgName")
        .arg(&contig)
        .arg("--ctgStart")
        .arg(start.to_string())
        .arg("--ctgEnd")
        .arg(end.to_string())
        .arg("--platform")
        .arg(&platform)
        .arg("--min_af")
        .arg("0.08")
        .arg("--indel_min_af")
        .arg("0.08")
        .arg("--extend_bed")
        .arg(&extended)
        .arg("--bed_fn")
        .arg(&bed)
        .arg("--unify_repre")
        .arg("--phasing_info_in_bam")
        .arg("--unify_repre_fn")
        .arg(&reference_rows)
        .output()
        .expect("failed to run the reference candidate discovery");
    assert!(
        status.status.success(),
        "the reference candidate discovery failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );

    let reference = parse_candidates(&fs::read_to_string(&reference_rows).unwrap());
    assert!(
        !reference.is_empty(),
        "the reference produced no candidates"
    );

    // Native candidate rows for the same region.
    let reference_fasta = Reference::open(&reference_path).unwrap();
    let mut config = CallConfig::new(
        bam,
        reference_path,
        PathBuf::new(),
        contig,
        if platform == "ilmn" {
            Platform::Ilmn
        } else if platform == "hifi" {
            Platform::Hifi
        } else {
            Platform::Ont
        },
    );
    config.start = Some(start);
    config.end = Some(end);
    config.bed = Some(bed);
    config.extend_bed = Some(extended);
    config.min_af = Some(0.08);
    let native = parse_candidates(&candidate_details_from_bam(&config, &reference_fasta).unwrap());

    let reference_only: Vec<_> = reference
        .keys()
        .filter(|position| !native.contains_key(position))
        .collect();
    let native_only: Vec<_> = native
        .keys()
        .filter(|position| !reference.contains_key(position))
        .collect();
    eprintln!(
        "candidate discovery: reference={} native={} shared={} reference-only={} native-only={}",
        reference.len(),
        native.len(),
        reference.len() - reference_only.len(),
        reference_only.len(),
        native_only.len()
    );
    eprintln!(
        "reference-only candidates: {:?}",
        &reference_only[..reference_only.len().min(8)]
    );
    eprintln!(
        "native-only candidates: {:?}",
        &native_only[..native_only.len().min(8)]
    );
    for position in reference_only.iter().take(4) {
        eprintln!(
            "reference {position}: {:?}",
            reference
                .get(position)
                .map(|(_, alleles)| alleles.keys().collect::<Vec<_>>())
        );
    }
    for position in native_only.iter().take(4) {
        eprintln!(
            "native {position}: {:?}",
            native
                .get(position)
                .map(|(_, alleles)| alleles.keys().collect::<Vec<_>>())
        );
    }

    // Every candidate that both implementations report must agree exactly in
    // depth and in the reads supporting each allele.
    for (position, (reference_depth, reference_alleles)) in &reference {
        let Some((native_depth, native_alleles)) = native.get(position) else {
            continue;
        };
        assert_eq!(
            reference_depth, native_depth,
            "candidate depth differs at {position}"
        );
        assert_eq!(
            reference_alleles, native_alleles,
            "candidate alleles differ at {position}"
        );
    }
}
