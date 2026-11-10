use std::{
    ffi::{OsStr, OsString},
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

const CHUNK_COUNT: usize = 8;

#[derive(Debug, Clone, PartialEq, Eq)]
struct UnifiedRecord {
    contig: String,
    position: u64,
    reference: String,
    alternate: String,
    genotype: String,
    genotype_quality: String,
    depth: String,
    allele_frequency: String,
    phase_set: Option<String>,
}

struct RealInputs {
    python: PathBuf,
    original: PathBuf,
    bam: PathBuf,
    truth: PathBuf,
    reference: PathBuf,
    bed: PathBuf,
    contig: String,
    platform: String,
}

fn run(command: &mut Command) -> Output {
    let display = format!("{command:?}");
    let output = command
        .output()
        .unwrap_or_else(|error| panic!("failed to run {display}: {error}"));
    assert!(
        output.status.success(),
        "{display} failed with {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn env_path(name: &str, default: &str) -> PathBuf {
    std::env::var_os(name)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default))
}

fn real_inputs(test_name: &str) -> Option<RealInputs> {
    let python = std::env::var_os("CLAIR3_ORIGINAL_PYTHON")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from("/Users/Aqua/Temp/Clair3/target/tmp/macos-validation-venv/bin/python")
        });
    let original = std::env::var_os("CLAIR3_ORIGINAL_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("/Users/Aqua/Temp/Clair3/target/tmp/original-main"));
    let inputs = RealInputs {
        python,
        original: original.join("clair3.py"),
        bam: env_path(
            "CLAIR3_UNIFICATION_BAM",
            "/Users/Aqua/Temp/Clair3/target/tmp/training-original-ru-inputs/chr21.ht.bam",
        ),
        truth: env_path(
            "CLAIR3_UNIFICATION_TRUTH",
            "/Users/Aqua/Temp/Clair3/target/tmp/training-original-ru-inputs/var_chr21",
        ),
        reference: env_path(
            "CLAIR3_UNIFICATION_REFERENCE",
            "/Users/Aqua/Temp/clairBenchmark/GRCh38_no_alt_analysis_set.fasta",
        ),
        bed: env_path(
            "CLAIR3_UNIFICATION_BED",
            "/Users/Aqua/Temp/Clair3/target/tmp/training-chr21-inputs/truth/region.bed",
        ),
        contig: std::env::var("CLAIR3_UNIFICATION_CONTIG").unwrap_or_else(|_| "chr21".to_string()),
        platform: std::env::var("CLAIR3_UNIFICATION_PLATFORM")
            .unwrap_or_else(|_| "hifi".to_string()),
    };
    let required = [
        ("validation Python", &inputs.python),
        ("original Clair3 entry point", &inputs.original),
        ("haplotagged chr21 BAM", &inputs.bam),
        ("original chr21 truth variants", &inputs.truth),
        ("GRCh38 reference FASTA", &inputs.reference),
        ("chr21 truth BED", &inputs.bed),
    ];
    let missing: Vec<_> = required
        .into_iter()
        .filter(|(_, path)| !path.is_file())
        .map(|(name, path)| format!("{name} ({})", path.display()))
        .collect();
    if !missing.is_empty() {
        eprintln!(
            "SKIP {test_name}: missing required input(s): {}",
            missing.join(", ")
        );
        return None;
    }
    Some(inputs)
}

fn unified_records(path: &Path) -> Vec<UnifiedRecord> {
    if !path.exists() {
        return Vec::new();
    }
    fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            assert!(fields.len() >= 10, "invalid unified VCF row: {line}");
            let format: Vec<_> = fields[8].split(':').collect();
            let sample: Vec<_> = fields[9].split(':').collect();
            let value = |name: &str| {
                let index = format
                    .iter()
                    .position(|field| *field == name)
                    .unwrap_or_else(|| panic!("missing FORMAT/{name} in row: {line}"));
                sample
                    .get(index)
                    .unwrap_or_else(|| {
                        panic!("missing sample value for FORMAT/{name} in row: {line}")
                    })
                    .to_string()
            };
            UnifiedRecord {
                contig: fields[0].to_string(),
                position: fields[1]
                    .parse()
                    .unwrap_or_else(|error| panic!("invalid VCF position in row {line}: {error}")),
                reference: fields[3].to_string(),
                alternate: fields[4].to_string(),
                genotype: value("GT"),
                genotype_quality: value("GQ"),
                depth: value("DP"),
                allele_frequency: value("AF"),
                phase_set: format.iter().position(|field| *field == "PS").map(|index| {
                    sample
                        .get(index)
                        .unwrap_or_else(|| {
                            panic!("missing sample value for FORMAT/PS in row: {line}")
                        })
                        .to_string()
                }),
            }
        })
        .collect()
}

fn unification_args(
    inputs: &RealInputs,
    extended_bed: &Path,
    output: &Path,
    chunk_id: usize,
    candidate_details: Option<&Path>,
    range: Option<(i64, i64)>,
) -> Vec<OsString> {
    let mut args = vec![
        OsString::from("UnifyRepresentation"),
        OsString::from("--platform"),
        OsString::from(inputs.platform.clone()),
        OsString::from("--bam_fn"),
        inputs.bam.as_os_str().to_owned(),
        OsString::from("--var_fn"),
        inputs.truth.as_os_str().to_owned(),
        OsString::from("--ref_fn"),
        inputs.reference.as_os_str().to_owned(),
        OsString::from("--bed_fn"),
        inputs.bed.as_os_str().to_owned(),
        OsString::from("--extend_bed"),
        extended_bed.as_os_str().to_owned(),
        OsString::from("--output_vcf_fn"),
        output.as_os_str().to_owned(),
        OsString::from("--min_af"),
        OsString::from("0.08"),
        OsString::from("--chunk_id"),
        OsString::from(chunk_id.to_string()),
        OsString::from("--chunk_num"),
        OsString::from(CHUNK_COUNT.to_string()),
        OsString::from("--ctgName"),
        OsString::from(inputs.contig.clone()),
    ];
    if let Some(candidate_details) = candidate_details {
        args.push(OsString::from("--candidate_details_fn"));
        args.push(candidate_details.as_os_str().to_owned());
    }
    if let Some((start, end)) = range {
        args.push(OsString::from("--ctgStart"));
        args.push(OsString::from(start.to_string()));
        args.push(OsString::from("--ctgEnd"));
        args.push(OsString::from(end.to_string()));
    }
    args.shrink_to_fit();
    args
}

fn run_original(inputs: &RealInputs, args: &[OsString]) {
    let mut command = Command::new(&inputs.python);
    command
        .arg(&inputs.original)
        .args(args)
        .arg("--pypy")
        .arg(&inputs.python)
        .env("PYTHONDONTWRITEBYTECODE", "1");
    run(&mut command);
}

fn run_native(args: &[OsString]) {
    run(Command::new(env!("CARGO_BIN_EXE_clair3")).args(args));
}

fn split_extend_bed(inputs: Option<&RealInputs>, bed: &Path, output: &Path) {
    let contig = inputs.map_or("chr21", |inputs| inputs.contig.as_str());
    let args: [&OsStr; 8] = [
        OsStr::new("SplitExtendBed"),
        OsStr::new("--bed_fn"),
        bed.as_os_str(),
        OsStr::new("--output_fn"),
        output.as_os_str(),
        OsStr::new("--ctgName"),
        OsStr::new(contig),
        OsStr::new("--expand_region_size"),
    ];
    if let Some(inputs) = inputs {
        let mut command = Command::new(&inputs.python);
        command
            .arg(&inputs.original)
            .args(args)
            .arg("33")
            .env("PYTHONDONTWRITEBYTECODE", "1");
        run(&mut command);
    } else {
        run(Command::new(env!("CARGO_BIN_EXE_clair3"))
            .args(args)
            .arg("33"));
    }
}

/// The region covered by the extended BED.
fn extended_bed_range(path: &Path) -> (i64, i64) {
    let mut start = i64::MAX;
    let mut end = i64::MIN;
    for line in fs::read_to_string(path).unwrap().lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 3 {
            start = start.min(fields[1].parse::<i64>().unwrap());
            end = end.max(fields[2].parse::<i64>().unwrap());
        }
    }
    assert!(start <= end, "the extended BED {path:?} is empty");
    (start, end)
}

/// The first and last candidate position of a candidate row file.
fn candidate_positions(path: &Path) -> (i64, i64) {
    let mut first = i64::MAX;
    let mut last = i64::MIN;
    for line in fs::read_to_string(path).unwrap().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let position: i64 = line
            .split_whitespace()
            .nth(1)
            .unwrap()
            .parse()
            .unwrap_or_else(|error| panic!("invalid candidate position in {line}: {error}"));
        first = first.min(position);
        last = last.max(position);
    }
    assert!(first <= last, "the candidate file {path:?} is empty");
    (first, last)
}

/// Generate the candidate rows with the reference implementation.
///
/// The reference applies an indel allele-frequency threshold inside its
/// candidate discovery, separately from `--min_af`, and the reference
/// unification command does not forward an indel threshold. The harness
/// therefore drives the reference candidate discovery directly with the same
/// thresholds the native implementation uses, so both sides work on the same
/// candidates.
fn original_candidate_rows(
    inputs: &RealInputs,
    extended_bed: &Path,
    output: &Path,
    start: i64,
    end: i64,
) {
    let mut command = Command::new(&inputs.python);
    command
        .arg(&inputs.original)
        .arg("CreateTensorFullAlignment")
        .arg("--bam_fn")
        .arg(&inputs.bam)
        .arg("--ref_fn")
        .arg(&inputs.reference)
        .arg("--ctgName")
        .arg(&inputs.contig)
        .arg("--ctgStart")
        .arg(start.to_string())
        .arg("--ctgEnd")
        .arg(end.to_string())
        .arg("--platform")
        .arg(&inputs.platform)
        .arg("--min_af")
        .arg("0.08")
        .arg("--indel_min_af")
        .arg("0.08")
        .arg("--extend_bed")
        .arg(extended_bed)
        .arg("--bed_fn")
        .arg(&inputs.bed)
        .arg("--unify_repre")
        .arg("--phasing_info_in_bam")
        .arg("--unify_repre_fn")
        .arg(output)
        .env("PYTHONDONTWRITEBYTECODE", "1");
    run(&mut command);
    assert!(
        output.exists() && fs::metadata(output).unwrap().len() > 0,
        "the reference produced no candidate rows"
    );
}

/// Compare the reference and native unification records.
///
/// The reference emits the alternate alleles of a multi-allelic record in
/// Python hash order, which is not reproducible between runs, while the native
/// implementation emits them in the candidate order. Every other field and the
/// allele set must match.
fn records_match(original: &[UnifiedRecord], native: &[UnifiedRecord]) -> Option<String> {
    let key = |record: &UnifiedRecord| (record.contig.clone(), record.position);
    let left: std::collections::BTreeMap<_, _> = original.iter().map(|r| (key(r), r)).collect();
    let right: std::collections::BTreeMap<_, _> = native.iter().map(|r| (key(r), r)).collect();
    if let Some(position) = left.keys().find(|position| !right.contains_key(position)) {
        return Some(format!("{position:?} is missing from the native output"));
    }
    if let Some(position) = right.keys().find(|position| !left.contains_key(position)) {
        return Some(format!("{position:?} is missing from the reference output"));
    }
    for ((_position, left), (_other, right)) in left.iter().zip(right.iter()) {
        for (name, left_value, right_value) in [
            ("reference allele", &left.reference, &right.reference),
            ("genotype", &left.genotype, &right.genotype),
            (
                "genotype quality",
                &left.genotype_quality,
                &right.genotype_quality,
            ),
            ("depth", &left.depth, &right.depth),
            (
                "allele frequency",
                &left.allele_frequency,
                &right.allele_frequency,
            ),
        ] {
            if left_value != right_value {
                return Some(format!(
                    "{name} differs: reference={left_value} native={right_value}"
                ));
            }
        }
        if left.phase_set != right.phase_set {
            return Some(format!(
                "phase set differs: reference={:?} native={:?}",
                left.phase_set, right.phase_set
            ));
        }
        let left_alleles: std::collections::BTreeSet<&str> = left.alternate.split(',').collect();
        let right_alleles: std::collections::BTreeSet<&str> = right.alternate.split(',').collect();
        if left_alleles != right_alleles {
            return Some(format!(
                "alternate alleles differ: reference={} native={}",
                left.alternate, right.alternate
            ));
        }
    }
    None
}

#[test]
fn split_extend_bed_keeps_final_interval() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("extended.bed");
    run(Command::new(env!("CARGO_BIN_EXE_clair3"))
        .arg("SplitExtendBed")
        .arg("--bed_fn")
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data/final_interval.bed"))
        .arg("--output_fn")
        .arg(&output)
        .arg("--expand_region_size")
        .arg("10")
        .arg("--ctgName")
        .arg("chrToy"));
    assert_eq!(
        fs::read_to_string(output).unwrap(),
        "chrToy 90 135\nchrToy 190 220"
    );
}

#[test]
#[ignore = "requires prepared original chr21 unification inputs"]
fn anchor_records_preserve_candidate_alleles_and_metadata() {
    let Some(inputs) = real_inputs("anchor_records_preserve_candidate_alleles_and_metadata") else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let extended = directory.path().join("original-extended.bed");
    split_extend_bed(Some(&inputs), &inputs.bed, &extended);
    let mut records = std::collections::BTreeMap::new();
    for chunk_id in 1..=2 {
        let output = directory.path().join(format!("rust-anchor-{chunk_id}.vcf"));
        run_native(&unification_args(
            &inputs, &extended, &output, chunk_id, None, None,
        ));
        for record in unified_records(&output) {
            records.insert(record.position, record);
        }
    }
    let expected = |position, reference: &str, alternate: &str, genotype: &str| UnifiedRecord {
        contig: "chr21".to_string(),
        position,
        reference: reference.to_string(),
        alternate: alternate.to_string(),
        genotype: genotype.to_string(),
        genotype_quality: "10".to_string(),
        depth: "10".to_string(),
        allele_frequency: "0.5000".to_string(),
        phase_set: None,
    };
    for expected in [
        expected(13048298, "TG", "T", "1/1"),
        expected(13093614, "CA", "CAA,C", "1/2"),
        expected(
            13245144,
            "GAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGA",
            "GAAGAAAGAAAGA,G",
            "1/2",
        ),
        expected(13251076, "C", "T", "1/1"),
    ] {
        assert_eq!(
            records.get(&expected.position),
            Some(&expected),
            "anchor chr21:{}",
            expected.position
        );
    }
    assert!(
        !records.contains_key(&13048299),
        "chr21:13048299 is a nonvariant in the matched original partition: {:?}",
        records.get(&13048299)
    );
}

#[test]
fn chunked_unification_matches_original() {
    let Some(inputs) = real_inputs("chunked_unification_matches_original") else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let original_extended = directory.path().join("original-extended.bed");
    let native_extended = directory.path().join("rust-extended.bed");
    split_extend_bed(Some(&inputs), &inputs.bed, &original_extended);
    split_extend_bed(None, &inputs.bed, &native_extended);
    let (start, end) = extended_bed_range(&original_extended);

    // Reference candidate rows with the same allele-frequency thresholds the
    // native implementation uses.
    let candidates = directory.path().join("original-candidates.txt");
    original_candidate_rows(&inputs, &original_extended, &candidates, start, end);

    // The reference and the native must agree on every record. The chunk
    // ranges are passed explicitly so both implementations cover exactly the
    // same region.
    let (first, last) = candidate_positions(&candidates);
    let size = (last - first) / CHUNK_COUNT as i64 + 1;
    let mut differences = Vec::new();
    let mut original_total = 0;
    let mut native_total = 0;
    for chunk_id in 1..=CHUNK_COUNT {
        let chunk_start = first + size * (chunk_id as i64 - 1);
        let chunk_end = chunk_start + size;
        let original_output = directory.path().join(format!("original-{chunk_id}.vcf"));
        let native_output = directory.path().join(format!("rust-{chunk_id}.vcf"));
        run_original(
            &inputs,
            &unification_args(
                &inputs,
                &original_extended,
                &original_output,
                chunk_id,
                Some(&candidates),
                Some((chunk_start, chunk_end)),
            ),
        );
        run_native(&unification_args(
            &inputs,
            &original_extended,
            &native_output,
            chunk_id,
            Some(&candidates),
            Some((chunk_start, chunk_end)),
        ));
        let original = unified_records(&original_output);
        let native = unified_records(&native_output);
        original_total += original.len();
        native_total += native.len();
        if let Some(difference) = records_match(&original, &native) {
            differences.push(format!("chunk {chunk_id}: {difference}"));
        }
    }
    assert!(
        original_total > 0 && native_total > 0,
        "empty unification comparison: reference records={original_total}, rust records={native_total}"
    );
    assert!(
        differences.is_empty(),
        "unification differed in {} of {CHUNK_COUNT} chunks; reference records={original_total}, rust records={native_total}\n{}",
        differences.len(),
        differences.join("\n")
    );
}
