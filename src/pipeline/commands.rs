use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufRead, BufWriter, Read as _, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use anyhow::{Context, Result, bail};
use clap::Args;

use crate::{
    config::Platform,
    io::{
        compression::open_text,
        reference::Reference,
        variant::{VcfRecord, vcf_header},
    },
    pipeline::{
        merge::{MergeOptions, merge_records, read_vcf, sort_records},
        metrics::{format_metrics, overall_metrics, overall_report},
        postprocess::{
            add_back_missing_variants, add_pair_end_depth, before_chrom_header, switch_zygosity,
        },
        run::{RunOptions, check_environment, compress_index_with_tools},
    },
};

#[derive(Args, Clone)]
#[command(rename_all = "snake_case")]
pub struct CompressionArgs {
    #[arg(long, default_value = "bgzip", hide = true)]
    bgzip: String,
    #[arg(long, default_value = "tabix", hide = true)]
    tabix: String,
    #[arg(long, default_value = "lz4", hide = true)]
    lz4: String,
}

fn uncompressed_output(path: &Path) -> PathBuf {
    if path.extension().and_then(|extension| extension.to_str()) == Some("gz") {
        path.with_extension("")
    } else {
        path.to_path_buf()
    }
}

fn write_vcf(path: &Path, header: &str, records: &[VcfRecord]) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let mut writer = BufWriter::new(fs::File::create(path)?);
    writer.write_all(header.as_bytes())?;
    for record in records {
        writeln!(writer, "{}", record.to_line())?;
    }
    writer.flush()?;
    Ok(())
}

fn write_compressed_vcf(
    output: &Path,
    header: &str,
    records: &[VcfRecord],
    tools: &CompressionArgs,
) -> Result<()> {
    let path = uncompressed_output(output);
    write_vcf(&path, header, records)?;
    compress_index_with_tools(&path, &tools.bgzip, &tools.tabix)
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct CheckEnvsArgs {
    #[arg(long = "bam_fn")]
    bam: PathBuf,
    #[arg(long = "ref_fn", default_value = "ref.fa")]
    reference: PathBuf,
    #[arg(long = "output_fn_prefix")]
    output: PathBuf,
    #[arg(long = "ctg_name", default_value = "EMPTY")]
    contig_names: String,
    #[arg(long = "bed_fn")]
    bed: Option<String>,
    #[arg(long = "vcf_fn")]
    known_vcf: Option<String>,
    #[arg(long, default_value_t = 5_000_000)]
    chunk_size: i64,
    #[arg(long, default_value_t = 0)]
    chunk_num: isize,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    include_all_ctgs: bool,
    #[arg(long, default_value_t = 16)]
    threads: usize,
    #[arg(long, default_value = "samtools")]
    samtools: String,
    #[arg(long, default_value = "whatshap")]
    whatshap: String,
    #[arg(long = "sampleName", alias = "sample_name", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long)]
    qual: Option<f64>,
    #[arg(long, default_value_t = 0.3)]
    var_pct_full: f64,
    #[arg(long, default_value_t = 0.3)]
    ref_pct_full: f64,
    #[arg(long, default_value_t = 0.08)]
    snp_min_af: f64,
    #[arg(long, default_value_t = 0.08)]
    indel_min_af: f64,
    #[arg(long, default_value_t = 0)]
    min_contig_size: u64,
    #[arg(long = "cmd_fn")]
    command_file: Option<PathBuf>,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    no_phasing_for_fa: bool,
    #[arg(long, hide = true)]
    python: Option<String>,
    #[arg(long, hide = true)]
    pypy: Option<String>,
    #[arg(long, hide = true)]
    parallel: Option<String>,
    #[command(flatten)]
    tools: CompressionArgs,
}

pub fn check_envs(args: CheckEnvsArgs) -> Result<()> {
    let _ = (args.python, args.pypy, args.parallel);
    let mut options = RunOptions::new(
        args.bam,
        args.reference,
        PathBuf::new(),
        args.output,
        Platform::Ont,
    );
    options.threads = args.threads;
    options.contig_names = (args.contig_names != "EMPTY").then_some(args.contig_names);
    options.bed = args.bed.filter(|path| path != "EMPTY").map(PathBuf::from);
    options.known_vcf = args
        .known_vcf
        .filter(|path| path != "EMPTY")
        .map(PathBuf::from);
    options.chunk_size = args.chunk_size;
    options.chunk_num = args.chunk_num;
    options.include_all_contigs = args.include_all_ctgs;
    options.samtools = args.samtools;
    options.whatshap = args.whatshap;
    options.bgzip = args.tools.bgzip;
    options.tabix = args.tools.tabix;
    options.sample_name = args.sample_name;
    options.qual = args.qual;
    options.var_pct_full = Some(args.var_pct_full);
    options.ref_pct_full = Some(args.ref_pct_full);
    options.snp_min_af = Some(args.snp_min_af);
    options.indel_min_af = Some(args.indel_min_af);
    options.min_contig_size = args.min_contig_size;
    options.no_phasing_for_fa = args.no_phasing_for_fa;
    options.command = args.command_file.map(fs::read_to_string).transpose()?;
    let contigs = check_environment(&options)?;
    println!("[INFO] Call variant in contigs: {}", contigs.join(" "));
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct SortVcfArgs {
    #[arg(long)]
    output_fn: PathBuf,
    #[arg(long)]
    input_dir: Option<PathBuf>,
    #[arg(long)]
    vcf_fn_prefix: Option<String>,
    #[arg(long, default_value = ".vcf")]
    vcf_fn_suffix: String,
    #[arg(long)]
    ref_fn: Option<PathBuf>,
    #[arg(long = "sampleName", alias = "sample_name", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long)]
    contigs_fn: Option<PathBuf>,
    #[arg(long = "cmd_fn")]
    command_file: Option<PathBuf>,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    pileup_only: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    print_ref_calls: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = true)]
    check_malformed_records: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_precise: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_sensitive: bool,
    #[arg(long, default_value = "2")]
    qual: Option<f64>,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    output_all_contigs_in_gvcf_header: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_gpu: bool,
    #[command(flatten)]
    tools: CompressionArgs,
}

fn contig_order(contigs: &[String], contig: &str) -> usize {
    let major = (1..=22)
        .map(|number| format!("chr{number}"))
        .chain(["chrX".to_string(), "chrY".to_string()])
        .chain((1..=22).map(|number| number.to_string()))
        .chain(["X".to_string(), "Y".to_string()])
        .position(|name| name == contig)
        .unwrap_or(usize::MAX);
    if major != usize::MAX {
        return major;
    }
    contigs
        .iter()
        .position(|name| name == contig)
        .map(|position| 48 + position)
        .unwrap_or(usize::MAX)
}

fn add_record(line: &str, contig: Option<&str>, records: &mut Vec<VcfRecord>) -> Result<()> {
    let record = VcfRecord::from_line(line)?;
    if contig.is_none_or(|contig| contig == record.contig) {
        records.push(record);
    }
    Ok(())
}

fn add_input_line(
    line: &str,
    contig: Option<&str>,
    sample_name: &str,
    repair: bool,
    header: &mut Vec<String>,
    records: &mut Vec<VcfRecord>,
) -> Result<()> {
    if line.starts_with('#') {
        if line.starts_with("##commandline=") {
            return Ok(());
        }
        if repair && line.starts_with("#CHROM") && line.split('\t').count() != 10 {
            let marker = if line.contains(sample_name) {
                sample_name
            } else {
                contig.context("malformed VCF header has no recognizable sample name")?
            };
            let split = line
                .find(marker)
                .map(|position| position + marker.len())
                .context("malformed VCF header has no sample boundary")?;
            let first = &line[..split];
            let second = &line[split..];
            if second.split('\t').count() > 10 {
                return Ok(());
            }
            if !header.iter().any(|row| row == first) {
                header.push(first.to_string());
            }
            if !second.is_empty() {
                add_record(second, contig, records)?;
            }
        } else if !header.iter().any(|row| row == line) {
            header.push(line.to_string());
        }
        return Ok(());
    }
    if line.is_empty() {
        return Ok(());
    }
    if repair {
        let columns = line.split('\t').count();
        if columns > 19 {
            eprintln!("[WARNING] Skipping malformed VCF record with {columns} columns");
            return Ok(());
        }
        if columns > 10
            && let Some(contig) = contig.or_else(|| line.split('\t').next())
            && let Some(split) = line[contig.len()..]
                .find(contig)
                .map(|position| position + contig.len())
        {
            add_record(&line[..split], Some(contig), records)?;
            add_record(&line[split..], Some(contig), records)?;
            return Ok(());
        }
    }
    add_record(line, contig, records)
}

fn read_sort_input(
    path: &Path,
    args: &SortVcfArgs,
    header: &mut Vec<String>,
    records: &mut Vec<VcfRecord>,
) -> Result<()> {
    let is_lz4 = path.is_file()
        && fs::File::open(path)
            .and_then(|mut file| {
                let mut magic = [0; 4];
                file.read_exact(&mut magic).map(|_| magic)
            })
            .is_ok_and(|magic| magic == [0x04, 0x22, 0x4d, 0x18]);
    let mut process = None;
    let reader: Box<dyn BufRead> = if is_lz4 {
        let mut child = Command::new(&args.tools.lz4)
            .arg("-d")
            .arg("-c")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .with_context(|| format!("failed to decompress LZ4 input {}", path.display()))?;
        let stdout = child
            .stdout
            .take()
            .context("LZ4 decompressor has no output")?;
        process = Some(child);
        Box::new(std::io::BufReader::new(stdout))
    } else {
        open_text(path).with_context(|| format!("failed to read {}", path.display()))?
    };
    let result = (|| {
        for line in reader.lines() {
            add_input_line(
                line?.trim_end_matches(['\r', '\n']).to_string().as_str(),
                None,
                &args.sample_name,
                args.check_malformed_records,
                header,
                records,
            )?;
        }
        Ok(())
    })();
    if let Some(mut process) = process {
        let status = process.wait()?;
        if !status.success() {
            bail!("LZ4 decompression failed for {}", path.display());
        }
    }
    result
}

fn read_contigs(path: Option<&Path>, records: &[VcfRecord]) -> Result<Vec<String>> {
    if let Some(path) = path {
        return Ok(fs::read_to_string(path)?
            .lines()
            .filter(|line| !line.is_empty())
            .map(str::to_string)
            .collect());
    }
    Ok(records
        .iter()
        .map(|record| record.contig.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect())
}

pub fn sort_vcf(args: SortVcfArgs) -> Result<()> {
    let _ = args.use_gpu;
    let mut header = Vec::new();
    let mut records = Vec::new();
    if let Some(directory) = &args.input_dir {
        let mut inputs = fs::read_dir(directory)
            .with_context(|| format!("failed to list input directory {}", directory.display()))?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?;
        inputs.retain(|path| {
            path.is_file()
                && path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| {
                        args.vcf_fn_prefix
                            .as_ref()
                            .is_none_or(|prefix| name.starts_with(prefix))
                            && name.ends_with(&args.vcf_fn_suffix)
                    })
        });
        inputs.sort();
        for input in inputs {
            read_sort_input(&input, &args, &mut header, &mut records)?;
        }
    } else {
        read_sort_input(Path::new("PIPE"), &args, &mut header, &mut records)?;
    }
    let contigs = read_contigs(args.contigs_fn.as_deref(), &records)?;
    if args.input_dir.is_some() {
        let allowed: BTreeSet<_> = contigs.iter().collect();
        records.retain(|record| allowed.contains(&record.contig));
    }
    let deduplicated: BTreeMap<_, _> = records
        .into_iter()
        .map(|record| ((record.contig.clone(), record.position), record))
        .collect();
    let mut records: Vec<_> = deduplicated.into_values().collect();
    records.sort_by_key(|record| {
        (
            contig_order(&contigs, &record.contig),
            record.contig.clone(),
            record.position,
        )
    });
    if args.pileup_only {
        records = merge_records(
            &records,
            &[],
            &MergeOptions {
                print_ref_calls: args.print_ref_calls,
                haploid_precise: args.haploid_precise,
                haploid_sensitive: args.haploid_sensitive,
                qual: args.qual,
            },
        );
    }
    if header.is_empty() {
        let reference_path = args
            .ref_fn
            .as_ref()
            .context("--ref_fn is required for an empty VCF")?;
        let reference = Reference::open(reference_path)?;
        let contigs = if contigs.is_empty() {
            reference.contigs()?
        } else {
            contigs.clone()
        };
        header = vcf_header(
            &reference,
            None,
            &args.sample_name,
            args.vcf_fn_suffix.contains("gvcf"),
            &contigs,
        )?
        .lines()
        .map(str::to_string)
        .collect();
    }
    if !header.iter().any(|line| line.starts_with("##cmdline="))
        && let Some(command_file) = &args.command_file
        && command_file.is_file()
    {
        let command = fs::read_to_string(command_file)?;
        let command = command.trim();
        if !command.is_empty() {
            let position = if header.len() >= 3 {
                3
            } else {
                header.len() - 1
            };
            header.insert(position, format!("##cmdline={command}"));
        }
    }
    if args.vcf_fn_suffix == ".gvcf" && !args.output_all_contigs_in_gvcf_header {
        let allowed: BTreeSet<_> = contigs.iter().map(String::as_str).collect();
        header.retain(|line| {
            line.strip_prefix("##contig=<ID=")
                .and_then(|value| value.split([',', '>']).next())
                .is_none_or(|contig| allowed.contains(contig))
        });
    }
    let header = header.join("\n") + "\n";
    write_vcf(&args.output_fn, &header, &records)?;
    if args.vcf_fn_suffix != ".tmp.gvcf" {
        compress_index_with_tools(&args.output_fn, &args.tools.bgzip, &args.tools.tabix)?;
    }
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct AddBackMissingArgs {
    #[arg(long = "vcf_fn")]
    known_vcf: PathBuf,
    #[arg(long = "clair3_input_vcf_fn")]
    called_vcf: PathBuf,
    #[arg(long)]
    output_fn: PathBuf,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = true)]
    switch_genotype: bool,
    #[command(flatten)]
    tools: CompressionArgs,
}

pub fn add_back_missing(args: AddBackMissingArgs) -> Result<()> {
    let known = read_vcf(&args.known_vcf)?;
    let called = read_vcf(&args.called_vcf)?;
    let mut records =
        add_back_missing_variants(&known.records, &called.records, args.switch_genotype)?;
    let contigs = read_contigs(None, &records)?;
    sort_records(&mut records, &contigs);
    println!(
        "[INFO] Total variants for genotyping: {}, total Clair3 variant output: {}, added {} variants into output VCF",
        known.records.len(),
        called.records.len(),
        known.records.len().saturating_sub(
            records
                .iter()
                .filter(|record| called
                    .records
                    .iter()
                    .any(|called| called.contig == record.contig
                        && called.position == record.position))
                .count()
        )
    );
    write_compressed_vcf(
        &args.output_fn,
        &(called.header.join("\n") + "\n"),
        &records,
        &args.tools,
    )
}

fn optional_tag(tag: Option<String>) -> Option<String> {
    tag.filter(|tag| !tag.eq_ignore_ascii_case("NONE"))
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct PairEndDepthArgs {
    #[arg(long = "bam_fn")]
    bam: PathBuf,
    #[arg(long = "ctg_name")]
    contig: Option<String>,
    #[arg(long = "clair3_vcf_input")]
    input_vcf: PathBuf,
    #[arg(long = "vcf_output")]
    output_vcf: PathBuf,
    #[arg(long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value = "PASS")]
    filter_tag: Option<String>,
    #[arg(long, default_value = "samtools", hide = true)]
    samtools: String,
    #[arg(long, default_value_t = 5)]
    min_mq: u8,
    #[arg(long, default_value_t = 0)]
    min_bq: u8,
    #[command(flatten)]
    tools: CompressionArgs,
}

pub fn pair_end_depth(args: PairEndDepthArgs) -> Result<()> {
    let _ = args.samtools;
    let input = read_vcf(&args.input_vcf)?;
    let records = add_pair_end_depth(
        &args.bam,
        &input.records,
        args.contig.as_deref(),
        optional_tag(args.filter_tag).as_deref(),
        args.threads,
        args.min_mq,
        args.min_bq,
    )?;
    let header = before_chrom_header(
        &input.header,
        &[
            "##FORMAT=<ID=PEAD,Number=1,Type=Integer,Description=\"Allelic depths for the alt alleles after removing same pair-end reads in the order listed\">",
        ],
    );
    write_compressed_vcf(&args.output_vcf, &header, &records, &args.tools)
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct SwitchZygosityArgs {
    #[arg(long = "bam_fn")]
    bam: PathBuf,
    #[arg(long = "ctg_name")]
    contig: Option<String>,
    #[arg(long = "clair3_vcf_input")]
    input_vcf: PathBuf,
    #[arg(long = "vcf_output")]
    output_vcf: PathBuf,
    #[arg(long = "sv_vcf_input")]
    sv_vcf: Option<PathBuf>,
    #[arg(long = "sv_bed_input")]
    sv_bed: Option<PathBuf>,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = true)]
    use_sv_qual: bool,
    #[arg(long, default_value_t = 8)]
    threads: usize,
    #[arg(long)]
    filter_tag: Option<String>,
    #[arg(long, default_value = "PASS")]
    sv_filter_tag: String,
    #[arg(long, default_value = "<DEL>")]
    sv_alt_tag: String,
    #[arg(long, default_value = "samtools", hide = true)]
    samtools: String,
    #[arg(long, default_value_t = 5)]
    min_mq: u8,
    #[arg(long, default_value_t = 0)]
    min_bq: u8,
    #[arg(long, default_value_t = 0.7)]
    max_af_for_zygosity_switching: f64,
    #[arg(long, default_value_t = 0.3)]
    min_del_proportion_for_long_sv: f64,
    #[arg(long, default_value_t = 16)]
    flanking_window_size: i64,
    #[command(flatten)]
    tools: CompressionArgs,
}

pub fn switch_zygosity_command(args: SwitchZygosityArgs) -> Result<()> {
    let _ = args.samtools;
    let input = read_vcf(&args.input_vcf)?;
    let before = input
        .records
        .iter()
        .filter(|record| {
            record
                .samples
                .first()
                .is_some_and(|sample| sample.starts_with("1/1"))
        })
        .count();
    let records = switch_zygosity(
        &args.bam,
        &input.records,
        args.contig.as_deref(),
        args.sv_vcf.as_deref(),
        args.sv_bed.as_deref(),
        args.use_sv_qual,
        optional_tag(args.filter_tag).as_deref(),
        &args.sv_filter_tag,
        &args.sv_alt_tag,
        args.max_af_for_zygosity_switching,
        args.min_del_proportion_for_long_sv,
        args.flanking_window_size,
        args.threads,
        args.min_mq,
        args.min_bq,
    )?;
    let switched = records
        .iter()
        .filter(|record| {
            record
                .samples
                .first()
                .is_some_and(|sample| sample.starts_with("0/1"))
                && record.info.contains("SVBASEDHET")
        })
        .count();
    let switched = if args.use_sv_qual && args.sv_bed.is_none() {
        switched
    } else {
        before.saturating_sub(
            records
                .iter()
                .filter(|record| {
                    record
                        .samples
                        .first()
                        .is_some_and(|sample| sample.starts_with("1/1"))
                })
                .count(),
        )
    };
    println!("[INFO] {switched} SNP variant switched from HOM to HET");
    let header = if args.use_sv_qual {
        before_chrom_header(
            &input.header,
            &[
                "##INFO=<ID=SVBASEDHET,Number=0,Type=Flag,Description=\"Short variant zygosity switched from HOM to HET given an overlapping HET SV\">",
                "##INFO=<ID=ORG_CLAIR3_SCORE,Number=1,Type=Float,Description=\"The original Clair3 quality score\">",
            ],
        )
    } else {
        input.header.join("\n") + "\n"
    };
    write_compressed_vcf(&args.output_vcf, &header, &records, &args.tools)
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct OverallMetricsArgs {
    #[arg(long = "happy_vcf_fn")]
    happy_vcf: PathBuf,
    #[arg(long = "log_happy")]
    happy_log: Option<PathBuf>,
    #[arg(long = "ctgName", alias = "ctg_name")]
    contig: Option<String>,
    #[arg(long)]
    output_fn: Option<PathBuf>,
}

pub fn overall_metrics_command(args: OverallMetricsArgs) -> Result<()> {
    let report = overall_report(
        &args.happy_vcf,
        args.contig.as_deref(),
        args.happy_log.as_deref(),
    )?;
    if let Some(output) = args.output_fn {
        if let Some(parent) = output.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)?;
        }
        fs::write(output, report)?;
    } else {
        print!("{report}");
    }
    Ok(())
}

pub fn metrics_text(path: &Path, contig: Option<&str>) -> Result<String> {
    Ok(format_metrics(&overall_metrics(path, contig)?))
}

pub fn ensure_output_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

pub fn validate_record_samples(records: &[VcfRecord]) -> Result<()> {
    if records.iter().any(|record| record.samples.is_empty()) {
        bail!("VCF record has no sample column");
    }
    Ok(())
}
