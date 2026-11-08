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

fn real_inputs() -> Option<RealInputs> {
    let inputs = RealInputs {
        python: PathBuf::from(
            "/Users/Aqua/Temp/Clair3/target/tmp/macos-validation-venv/bin/python",
        ),
        original: PathBuf::from("/Users/Aqua/Temp/Clair3/target/tmp/original-main/clair3.py"),
        bam: PathBuf::from(
            "/Users/Aqua/Temp/Clair3/target/tmp/training-original-ru-inputs/chr21.ht.bam",
        ),
        truth: PathBuf::from(
            "/Users/Aqua/Temp/Clair3/target/tmp/training-original-ru-inputs/var_chr21",
        ),
        reference: PathBuf::from(
            "/Users/Aqua/Temp/clairBenchmark/GRCh38_no_alt_analysis_set.fasta",
        ),
        bed: PathBuf::from(
            "/Users/Aqua/Temp/Clair3/target/tmp/training-chr21-inputs/truth/region.bed",
        ),
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
            "SKIP chunked_unification_matches_original: missing required input(s): {}",
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
) -> Vec<OsString> {
    let mut args = vec![
        OsString::from("UnifyRepresentation"),
        OsString::from("--platform"),
        OsString::from("hifi"),
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
        OsString::from("chr21"),
    ];
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
    let args: [&OsStr; 8] = [
        OsStr::new("SplitExtendBed"),
        OsStr::new("--bed_fn"),
        bed.as_os_str(),
        OsStr::new("--output_fn"),
        output.as_os_str(),
        OsStr::new("--ctgName"),
        OsStr::new("chr21"),
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

fn chunk_difference(
    chunk_id: usize,
    original: &[UnifiedRecord],
    native: &[UnifiedRecord],
) -> String {
    let mut positions: std::collections::BTreeMap<
        (String, u64),
        (Vec<&UnifiedRecord>, Vec<&UnifiedRecord>),
    > = std::collections::BTreeMap::new();
    for record in original {
        positions
            .entry((record.contig.clone(), record.position))
            .or_default()
            .0
            .push(record);
    }
    for record in native {
        positions
            .entry((record.contig.clone(), record.position))
            .or_default()
            .1
            .push(record);
    }
    let difference_count = positions
        .iter()
        .filter(|(_, (left, right))| left != right)
        .count();
    let mut differences = Vec::new();
    for ((contig, position), (left, right)) in
        positions.iter().filter(|(_, (left, right))| left != right)
    {
        differences.push(format!(
            "  {contig}:{position}: original={left:?}\n             rust={right:?}"
        ));
        if differences.len() == 4 {
            break;
        }
    }
    if differences.is_empty() {
        differences.push("  record order differs; per-position records are identical".to_string());
    }
    format!(
        "chunk {chunk_id}: original records={}, rust records={}, differing positions={difference_count}\n{}",
        original.len(),
        native.len(),
        differences.join("\n")
    )
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
    let Some(inputs) = real_inputs() else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let extended = directory.path().join("original-extended.bed");
    split_extend_bed(Some(&inputs), &inputs.bed, &extended);
    let mut records = std::collections::BTreeMap::new();
    for chunk_id in 1..=2 {
        let output = directory.path().join(format!("rust-anchor-{chunk_id}.vcf"));
        run_native(&unification_args(&inputs, &extended, &output, chunk_id));
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
#[ignore = "requires prepared original chr21 unification inputs"]
fn chunked_unification_matches_original() {
    let Some(inputs) = real_inputs() else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let original_extended = directory.path().join("original-extended.bed");
    let native_extended = directory.path().join("rust-extended.bed");
    split_extend_bed(Some(&inputs), &inputs.bed, &original_extended);
    split_extend_bed(None, &inputs.bed, &native_extended);

    let mut differences = Vec::new();
    let mut original_total = 0;
    let mut native_total = 0;
    for chunk_id in 1..=CHUNK_COUNT {
        let original_output = directory.path().join(format!("original-{chunk_id}.vcf"));
        let native_output = directory.path().join(format!("rust-{chunk_id}.vcf"));
        run_original(
            &inputs,
            &unification_args(&inputs, &original_extended, &original_output, chunk_id),
        );
        run_native(&unification_args(
            &inputs,
            &original_extended,
            &native_output,
            chunk_id,
        ));
        let original = unified_records(&original_output);
        let native = unified_records(&native_output);
        original_total += original.len();
        native_total += native.len();
        if original != native {
            differences.push(chunk_difference(chunk_id, &original, &native));
        }
    }
    assert!(
        original_total > 0 && native_total > 0,
        "empty unification comparison: original records={original_total}, rust records={native_total}"
    );
    assert!(
        differences.is_empty(),
        "unification differed in {} of {CHUNK_COUNT} chunks; original records={original_total}, rust records={native_total}\n{}",
        differences.len(),
        differences.join("\n")
    );
}
