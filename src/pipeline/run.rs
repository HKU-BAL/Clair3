use std::{
    collections::{BTreeMap, BTreeSet},
    ffi::OsStr,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
        mpsc,
    },
};

use anyhow::{Context, Result, bail};
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use tch::Device;

use crate::{
    call::caller::write_vcf,
    config::Platform,
    io::{
        interval::IntervalMap,
        reference::Reference,
        variant::{VcfRecord, vcf_header},
    },
    model::{Model, ModelKind},
    pipeline::{
        call::{
            CallConfig, CallOutput, Inference, InferenceRequest, call_full_alignment, call_pileup,
            checkpoint_path, inference_job,
        },
        gvcf::{GvcfOptions, merge_nonvariants, write_gvcf},
        merge::{MergeOptions, merge_illumina_records, merge_records, read_vcf, sort_records},
        phasing::{parallel_phase_command, phase_regions_from_files, scan_phase_regions},
        select::{SelectOptions, full_quality, phase_quality, select_candidates, select_het_snp},
    },
};

#[derive(Debug, Clone)]
pub struct RunOptions {
    pub bam: PathBuf,
    pub reference: PathBuf,
    pub model_path: PathBuf,
    pub output: PathBuf,
    pub platform: Platform,
    pub threads: usize,
    pub bed: Option<PathBuf>,
    pub known_vcf: Option<PathBuf>,
    pub contig_names: Option<String>,
    pub sample_name: String,
    pub chunk_num: isize,
    pub chunk_size: i64,
    pub qual: Option<f64>,
    pub samtools: String,
    pub whatshap: String,
    pub longphase: String,
    pub bgzip: String,
    pub tabix: String,
    pub var_pct_full: Option<f64>,
    pub ref_pct_full: Option<f64>,
    pub var_pct_phasing: Option<f64>,
    pub snp_min_af: Option<f64>,
    pub indel_min_af: Option<f64>,
    pub min_mq: u8,
    pub min_bq: u8,
    pub min_coverage: i64,
    pub min_contig_size: u64,
    pub pileup_model_prefix: String,
    pub full_alignment_model_prefix: String,
    pub pileup_only: bool,
    pub fast_mode: bool,
    pub call_snp_only: bool,
    pub call_head_tail: bool,
    pub output_all_contigs_in_gvcf_header: bool,
    pub print_ref_calls: bool,
    pub gvcf: bool,
    pub gq_bin_size: i32,
    pub base_error: f64,
    pub bp_resolution: bool,
    pub haploid_precise: bool,
    pub haploid_sensitive: bool,
    pub include_all_contigs: bool,
    pub no_phasing_for_fa: bool,
    pub remove_intermediate_dir: bool,
    pub enable_long_indel: bool,
    pub keep_iupac: bool,
    pub use_gpu: bool,
    pub device: String,
    pub use_longphase_for_intermediate_phasing: bool,
    pub use_whatshap_for_final_phasing: bool,
    pub use_longphase_for_final_phasing: bool,
    pub use_whatshap_for_final_haplotagging: bool,
    pub enable_dwell_time: bool,
    pub need_realignment: bool,
    pub command: Option<String>,
}

impl RunOptions {
    pub fn new(
        bam: PathBuf,
        reference: PathBuf,
        model_path: PathBuf,
        output: PathBuf,
        platform: Platform,
    ) -> Self {
        Self {
            bam,
            reference,
            model_path,
            output,
            platform,
            threads: 1,
            bed: None,
            known_vcf: None,
            contig_names: None,
            sample_name: "SAMPLE".to_string(),
            chunk_num: 0,
            chunk_size: 5_000_000,
            qual: Some(2.0),
            samtools: "samtools".to_string(),
            whatshap: "whatshap".to_string(),
            longphase: "longphase".to_string(),
            bgzip: "bgzip".to_string(),
            tabix: "tabix".to_string(),
            var_pct_full: None,
            ref_pct_full: None,
            var_pct_phasing: None,
            snp_min_af: None,
            indel_min_af: None,
            min_mq: 5,
            min_bq: 0,
            min_coverage: 2,
            min_contig_size: 0,
            pileup_model_prefix: "pileup".to_string(),
            full_alignment_model_prefix: "full_alignment".to_string(),
            pileup_only: false,
            fast_mode: false,
            call_snp_only: false,
            call_head_tail: false,
            output_all_contigs_in_gvcf_header: false,
            print_ref_calls: false,
            gvcf: false,
            gq_bin_size: 5,
            base_error: 0.001,
            bp_resolution: false,
            haploid_precise: false,
            haploid_sensitive: false,
            include_all_contigs: false,
            no_phasing_for_fa: false,
            remove_intermediate_dir: false,
            enable_long_indel: false,
            keep_iupac: false,
            use_gpu: false,
            device: "EMPTY".to_string(),
            use_longphase_for_intermediate_phasing: false,
            use_whatshap_for_final_phasing: false,
            use_longphase_for_final_phasing: false,
            use_whatshap_for_final_haplotagging: false,
            enable_dwell_time: false,
            need_realignment: true,
            command: None,
        }
    }
}

fn execute(command: &mut Command) -> Result<()> {
    let description = format!("{command:?}");
    let status = command
        .status()
        .with_context(|| format!("failed to execute {description}"))?;
    if !status.success() {
        bail!("external command failed with {status}: {description}");
    }
    Ok(())
}

pub fn validate_tool_version(executable: &str, name: &str, minimum: &[u64]) -> Result<()> {
    let output = Command::new(executable)
        .arg("--version")
        .output()
        .with_context(|| format!("failed to execute {name} at {executable}"))?;
    if !output.status.success() {
        bail!("{name} --version failed with {}", output.status);
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let first_line = stdout
        .lines()
        .next()
        .with_context(|| format!("{name} returned no version information"))?;
    let version_text = first_line
        .split_whitespace()
        .find(|field| field.starts_with(|character: char| character.is_ascii_digit()))
        .with_context(|| format!("could not parse {name} version from {first_line:?}"))?;
    let mut version = version_text
        .split('.')
        .map(|part| {
            part.parse::<u64>()
                .with_context(|| format!("could not parse {name} version {version_text:?}"))
        })
        .collect::<Result<Vec<_>>>()?;
    version.resize(minimum.len(), 0);
    if version.as_slice() < minimum {
        let required = minimum
            .iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join(".");
        bail!("{name} {version_text} is too old; version {required} or later is required");
    }
    Ok(())
}

fn tool_available(executable: &str) -> bool {
    Command::new(executable)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .output()
        .is_ok_and(|output| output.status.success())
}

fn validate_tools(options: &RunOptions) -> Result<()> {
    validate_tool_version(&options.samtools, "samtools", &[1, 15, 1])?;
    if options.pileup_only {
        return Ok(());
    }
    if !options.no_phasing_for_fa && options.platform != Platform::Ilmn {
        if options.use_longphase_for_intermediate_phasing {
            if !tool_available(&options.longphase) {
                bail!("LongPhase not found at {}", options.longphase);
            }
        } else if tool_available(&options.whatshap) {
            validate_tool_version(&options.whatshap, "whatshap", &[1, 7])?;
        } else if !tool_available(&options.longphase) {
            bail!(
                "phasing requires WhatsHap >= 1.7 or LongPhase; install one or use --no_phasing_for_fa"
            );
        }
    }
    if (options.platform == Platform::Ilmn && !options.no_phasing_for_fa)
        || options.use_whatshap_for_final_phasing
        || options.use_whatshap_for_final_haplotagging
    {
        validate_tool_version(&options.whatshap, "whatshap", &[1, 7])?;
    }
    if options.use_longphase_for_final_phasing && !tool_available(&options.longphase) {
        bail!("LongPhase not found at {}", options.longphase);
    }
    Ok(())
}

pub(crate) fn compress_index_with_tools(path: &Path, bgzip: &str, tabix: &str) -> Result<()> {
    execute(Command::new(bgzip).arg("-f").arg(path))?;
    let compressed = path.with_extension(format!(
        "{}.gz",
        path.extension().and_then(OsStr::to_str).unwrap_or("vcf")
    ));
    execute(
        Command::new(tabix)
            .arg("-f")
            .arg("-p")
            .arg("vcf")
            .arg(&compressed),
    )?;
    Ok(())
}

pub(crate) fn compress_index(path: &Path, options: &RunOptions) -> Result<()> {
    compress_index_with_tools(path, &options.bgzip, &options.tabix)
}

pub(crate) fn write_output_vcf(
    path: &Path,
    reference: &Reference,
    records: &[VcfRecord],
    contigs: &[String],
    options: &RunOptions,
    gvcf: bool,
) -> Result<()> {
    let mut header = vcf_header(
        reference,
        options.command.as_deref(),
        &options.sample_name,
        gvcf,
        &reference.contigs()?,
    )?;
    if gvcf && !options.output_all_contigs_in_gvcf_header {
        header = header
            .lines()
            .filter(|line| {
                !line.starts_with("##contig=")
                    || contigs
                        .iter()
                        .any(|contig| line.starts_with(&format!("##contig=<ID={contig},")))
            })
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut writer = BufWriter::new(fs::File::create(path)?);
    writer.write_all(header.as_bytes())?;
    for record in records {
        writeln!(writer, "{}", record.to_line())?;
    }
    writer.flush()?;
    compress_index(path, options)
}

fn major_contig_order(contig: &str) -> usize {
    let numeric = contig.strip_prefix("chr").unwrap_or(contig);
    let suffix = if let Ok(number) = numeric.parse::<usize>() {
        number
    } else if numeric == "X" {
        23
    } else if numeric == "Y" {
        24
    } else {
        return usize::MAX;
    };
    if contig.starts_with("chr") {
        suffix - 1
    } else {
        suffix + 23
    }
}

fn prepare_contigs(
    reference: &Reference,
    bam_path: &Path,
    options: &RunOptions,
) -> Result<Vec<String>> {
    let reference_contigs = reference.contigs()?;
    let bed_contigs: BTreeSet<String> = options
        .bed
        .as_ref()
        .map(|path| {
            IntervalMap::from_bed(path, None).map(|map| map.contigs().map(str::to_string).collect())
        })
        .transpose()?
        .unwrap_or_default();
    let known_contigs: BTreeSet<String> = options
        .known_vcf
        .as_ref()
        .map(|path| {
            read_vcf(path).map(|file| {
                file.records
                    .into_iter()
                    .map(|record| record.contig)
                    .collect()
            })
        })
        .transpose()?
        .unwrap_or_default();
    let requested: Option<BTreeSet<String>> = options.contig_names.as_ref().map(|names| {
        names
            .split(',')
            .filter(|name| !name.is_empty() && *name != "EMPTY")
            .map(str::to_string)
            .collect()
    });
    let mut reader = IndexedReader::from_path(bam_path)?;
    reader.set_reference(reference.path())?;
    let mut contigs = Vec::new();
    for contig in reference_contigs {
        if let Some(requested) = &requested
            && !requested.contains(&contig)
        {
            continue;
        }
        if options.bed.is_some() && !bed_contigs.contains(&contig)
            || options.known_vcf.is_some() && !known_contigs.contains(&contig)
        {
            continue;
        }
        if !options.include_all_contigs
            && requested.is_none()
            && options.bed.is_none()
            && options.known_vcf.is_none()
            && major_contig_order(&contig) == usize::MAX
        {
            continue;
        }
        if reference.contig_len(&contig) < options.min_contig_size {
            continue;
        }
        if reader.fetch(contig.as_bytes()).is_err() {
            continue;
        }
        if reader.records().next().transpose()?.is_none() {
            continue;
        }
        contigs.push(contig);
    }
    contigs.sort_by_key(|contig| (major_contig_order(contig), contig.clone()));
    Ok(contigs)
}

fn write_extended_regions(
    path: &Path,
    bed: Option<&PathBuf>,
    known_vcf: Option<&PathBuf>,
    contigs: &[String],
) -> Result<()> {
    let mut rows: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
    if let Some(path) = bed {
        for contig in contigs {
            let regions = IntervalMap::from_bed(path, Some(contig))?;
            for region in regions.regions(contig) {
                rows.entry(contig.clone())
                    .or_default()
                    .push(((region.start - 33).max(0), region.end + 33));
            }
        }
    } else if let Some(path) = known_vcf {
        for record in read_vcf(path)?.records {
            if record.position <= 33 {
                continue;
            }
            rows.entry(record.contig)
                .or_default()
                .push((record.position - 1 - 33, record.position + 33));
        }
    }
    fs::create_dir_all(path)?;
    for (contig, intervals) in rows {
        let content = intervals
            .into_iter()
            .map(|(start, end)| format!("{contig}\t{}\t{end}", start.max(0)))
            .collect::<Vec<_>>()
            .join("\n");
        fs::write(path.join(contig), content)?;
    }
    Ok(())
}

fn call_config(
    options: &RunOptions,
    checkpoint: &Path,
    bam: &Path,
    contig: &str,
    device: Device,
    var_pct: f64,
) -> CallConfig {
    let mut config = CallConfig::new(
        bam.to_path_buf(),
        options.reference.clone(),
        checkpoint.to_path_buf(),
        contig.to_string(),
        options.platform,
    );
    config.sample_name = options.sample_name.clone();
    config.bed = options.bed.clone();
    config.known_vcf = options.known_vcf.clone();
    config.min_mq = options.min_mq;
    config.min_bq = options.min_bq;
    config.min_coverage = options.min_coverage;
    config.snp_min_af = options.snp_min_af.unwrap_or(0.08);
    config.indel_min_af = options
        .indel_min_af
        .unwrap_or(options.platform.min_indel_af());
    config.fast_mode = options.fast_mode;
    config.call_snp_only = options.call_snp_only;
    config.call_head_tail = options.call_head_tail;
    config.enable_long_indel = options.enable_long_indel;
    config.enable_dwell_time = options.enable_dwell_time;
    config.haploid_precise = options.haploid_precise;
    config.haploid_sensitive = options.haploid_sensitive;
    config.no_phasing_for_fa = options.no_phasing_for_fa;
    config.need_realignment = options.need_realignment;
    config.gvcf = options.gvcf;
    config.gvcf_options = GvcfOptions {
        base_error: options.base_error,
        gq_bin_size: options.gq_bin_size,
        bp_resolution: options.bp_resolution,
    };
    config.device = device;
    config.quality_for_pass = options.qual;
    config.command = options.command.clone();
    config.keep_iupac = options.keep_iupac;
    config.min_af = (var_pct < 0.0).then_some(0.08);
    config
}

fn execute_call_jobs(jobs: &[CallConfig], threads: usize, pileup: bool) -> Result<Vec<CallOutput>> {
    if jobs.is_empty() {
        return Ok(Vec::new());
    }
    let count = jobs.len();
    let jobs = Arc::new(jobs);
    let next = Arc::new(AtomicUsize::new(0));
    let (sender, receiver) = mpsc::channel::<(usize, Result<CallOutput>)>();
    let workers = threads.max(1).min(count);
    let create_model = |config: &CallConfig| {
        Model::new(
            if pileup {
                ModelKind::Pileup
            } else {
                ModelKind::FullAlignment
            },
            config.device,
            Some(if pileup {
                18
            } else if config.enable_dwell_time {
                9
            } else {
                8
            }),
            !pileup && config.add_indel_length,
            Some(&config.checkpoint),
        )
    };
    let first = jobs.first().context("calling worker received no jobs")?;
    let (inference_sender, inference_receiver) =
        matches!(first.device, Device::Mps | Device::Cuda(_))
            .then(mpsc::channel::<InferenceRequest>)
            .unzip();
    std::thread::scope(|scope| {
        let inference_workers = inference_receiver
            .map(|receiver| {
                let count = std::env::var("CLAIR3_INFERENCE_WORKERS")
                    .ok()
                    .and_then(|value| value.parse::<usize>().ok())
                    .unwrap_or(1)
                    .max(1);
                let receiver = std::sync::Arc::new(std::sync::Mutex::new(receiver));
                (0..count)
                    .map(|_| {
                        let receiver = std::sync::Arc::clone(&receiver);
                        let config = first.clone();
                        scope.spawn(move || -> Result<()> {
                            let model = create_model(&config)?;
                            loop {
                                let request = receiver
                                    .lock()
                                    .map_err(|_| anyhow::anyhow!("inference queue lock poisoned"))?
                                    .recv();
                                let Ok((job, response)) = request else {
                                    break;
                                };
                                response.send(inference_job(job, &model)).map_err(|_| {
                                    anyhow::anyhow!("inference response channel closed")
                                })?;
                            }
                            Ok(())
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let mut handles = Vec::new();
        for _ in 0..workers {
            let jobs = Arc::clone(&jobs);
            let next = Arc::clone(&next);
            let sender = sender.clone();
            let inference_sender = inference_sender.clone();
            handles.push(scope.spawn(move || {
                let worker = (|| -> Result<()> {
                    let first = jobs.first().context("calling worker received no jobs")?;
                    let reference = Reference::open(&first.reference)?;
                    let model = inference_sender
                        .is_none()
                        .then(|| create_model(first))
                        .transpose()?
                        .map(std::sync::Mutex::new);
                    let inference = if let Some(sender) = &inference_sender {
                        Inference::Remote(sender)
                    } else {
                        Inference::Local(model.as_ref().unwrap())
                    };
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= count {
                            break;
                        }
                        let result = if pileup {
                            call_pileup(&jobs[index], &reference, &inference)
                        } else {
                            call_full_alignment(&jobs[index], &reference, &inference)
                        };
                        if sender.send((index, result)).is_err() {
                            break;
                        }
                    }
                    Ok(())
                })();
                if let Err(error) = worker {
                    let _ = sender.send((usize::MAX, Err(error)));
                }
            }));
        }
        drop(sender);
        let mut outputs = BTreeMap::new();
        let mut errors = Vec::new();
        for (index, result) in receiver {
            match result {
                Ok(output) => {
                    outputs.insert(index, output);
                }
                Err(error) => errors.push(error),
            }
        }
        for handle in handles {
            if handle.join().is_err() {
                bail!("a calling worker panicked");
            }
        }
        drop(inference_sender);
        for worker in inference_workers {
            worker
                .join()
                .map_err(|_| anyhow::anyhow!("inference worker panicked"))??;
        }
        if let Some(error) = errors.into_iter().next() {
            return Err(error);
        }
        if outputs.len() != count {
            bail!("calling workers returned {} of {count} jobs", outputs.len());
        }
        Ok(outputs.into_values().collect())
    })
}

fn link_bam(source: &Path, destination: &Path) -> Result<()> {
    let source = source.canonicalize()?;
    if destination.exists() || destination.is_symlink() {
        fs::remove_file(destination)?;
    }
    #[cfg(unix)]
    std::os::unix::fs::symlink(&source, destination)?;
    #[cfg(not(unix))]
    fs::copy(&source, destination)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn phase_bams(
    options: &RunOptions,
    reference: &Reference,
    contigs: &[String],
    pileup: &[VcfRecord],
    pileup_header: &[String],
    phase_vcf: &Path,
    phase_bam: &Path,
    phasing_proportion: f64,
    var_pct_full: f64,
) -> Result<()> {
    fs::create_dir_all(phase_vcf)?;
    fs::create_dir_all(phase_bam)?;
    let link_input_bam = |contig: &str| -> Result<()> {
        let bam = phase_bam.join(format!("{contig}.bam"));
        link_bam(&options.bam, &bam)?;
        for index in [
            options.bam.with_extension("bam.bai"),
            options.bam.with_extension("bai"),
            options.bam.with_extension("bam.csi"),
            options.bam.with_extension("csi"),
            options.bam.with_extension("bam.crai"),
            options.bam.with_extension("crai"),
        ] {
            if index.is_file() {
                link_bam(&index, &PathBuf::from(format!("{}.bai", bam.display())))?;
                break;
            }
        }
        Ok(())
    };
    if options.no_phasing_for_fa {
        for contig in contigs {
            link_input_bam(contig)?;
        }
        return Ok(());
    }
    let threshold = phase_quality(pileup, Some(phasing_proportion), var_pct_full);
    let use_longphase = options.platform != Platform::Ilmn
        && (options.use_longphase_for_intermediate_phasing
            || (!tool_available(&options.whatshap) && tool_available(&options.longphase)));
    fs::write(phase_vcf.join("phase_qual"), threshold.to_string())?;
    for contig in contigs {
        let selected = select_het_snp(pileup, contig, var_pct_full, Some(threshold));
        let path = phase_vcf.join(format!("{contig}.vcf"));
        let mut writer = BufWriter::new(fs::File::create(&path)?);
        writer.write_all(pileup_header.join("\n").as_bytes())?;
        writer.write_all(b"\n")?;
        for record in &selected {
            writeln!(writer, "{}", record.to_line())?;
        }
        writer.flush()?;
        let phased = phase_vcf.join(format!("phased_{contig}.vcf.gz"));
        if selected.is_empty() {
            let unphased = phase_vcf.join(format!("phased_{contig}.vcf"));
            fs::copy(&path, &unphased)?;
            compress_index(&unphased, options)?;
            link_input_bam(contig)?;
            continue;
        }
        if use_longphase {
            let prefix = phase_vcf.join(format!("phased_{contig}"));
            execute(
                Command::new(&options.longphase)
                    .arg("phase")
                    .arg("-s")
                    .arg(&path)
                    .arg("-b")
                    .arg(&options.bam)
                    .arg("-r")
                    .arg(reference.path())
                    .arg("-t")
                    .arg(options.threads.to_string())
                    .arg("-o")
                    .arg(&prefix)
                    .arg(if options.platform == Platform::Ont {
                        "--ont"
                    } else {
                        "--pb"
                    }),
            )?;
            execute(
                Command::new(&options.bgzip)
                    .arg("-f")
                    .arg(prefix.with_extension("vcf")),
            )?;
        } else {
            let mut phaser = Command::new(&options.whatshap);
            if options.platform != Platform::Ilmn {
                let positions: Vec<_> = selected
                    .iter()
                    .map(VcfRecord::position_zero_based)
                    .collect();
                let span_file = phase_vcf.join(format!("{contig}.spans"));
                let regions = if span_file.is_file() {
                    phase_regions_from_files(&positions, &[span_file])?
                } else {
                    scan_phase_regions(&options.bam, contig, &positions)?
                };
                let region_path = phase_vcf.join(format!("{contig}.phase-regions"));
                let mut region_writer = BufWriter::new(fs::File::create(&region_path)?);
                for region in regions {
                    writeln!(
                        region_writer,
                        "{}\t{}\t{}\t{}",
                        region.start, region.end, region.first_variant, region.end_variant
                    )?;
                }
                region_writer.flush()?;
                phaser =
                    parallel_phase_command(&options.whatshap, options.threads, Some(&region_path));
            }
            execute(
                phaser
                    .arg("phase")
                    .arg("--output")
                    .arg(&phased)
                    .arg("--reference")
                    .arg(reference.path())
                    .arg("--chromosome")
                    .arg(contig)
                    .arg("--distrust-genotypes")
                    .arg("--ignore-read-groups")
                    .arg(&path)
                    .arg(&options.bam),
            )?;
        }
        execute(
            Command::new(&options.tabix)
                .arg("-f")
                .arg("-p")
                .arg("vcf")
                .arg(&phased),
        )?;
        if options.platform != Platform::Ilmn {
            continue;
        }
        let bam = phase_bam.join(format!("{contig}.bam"));
        let mut haplotag = Command::new(&options.whatshap);
        haplotag
            .arg("haplotag")
            .arg("--output")
            .arg(&bam)
            .arg("--reference")
            .arg(reference.path())
            .arg("--ignore-read-groups");
        if !contig.contains(':') {
            haplotag.arg("--regions").arg(contig);
        }
        execute(haplotag.arg(&phased).arg(&options.bam))?;
        execute(
            Command::new(&options.samtools)
                .arg("index")
                .arg("-@")
                .arg(options.threads.to_string())
                .arg(&bam),
        )?;
    }
    Ok(())
}

fn add_back_missing(known: &Path, records: &mut Vec<VcfRecord>) -> Result<()> {
    let existing: BTreeSet<(String, i64)> = records
        .iter()
        .map(|record| (record.contig.clone(), record.position))
        .collect();
    for mut record in read_vcf(known)?.records {
        if existing.contains(&(record.contig.clone(), record.position)) {
            continue;
        }
        record.reference = record.reference.chars().next().unwrap_or('.').to_string();
        record.alternate = ".".to_string();
        record.quality = ".".to_string();
        record.filter = ".".to_string();
        record.info = ".".to_string();
        record.format = Some("GT".to_string());
        record.samples = vec!["./.".to_string()];
        records.push(record);
    }
    Ok(())
}

pub fn validate_dwell_tags(bam: &Path, reference: &Path) -> Result<()> {
    let mut reader = IndexedReader::from_path(bam)?;
    reader.set_reference(reference)?;
    reader.fetch(rust_htslib::bam::FetchDefinition::All)?;
    let mut checked = 0;
    let mut invalid = false;
    for record in reader.records().take(50) {
        let record = record?;
        checked += 1;
        match record.aux(b"mv") {
            Ok(rust_htslib::bam::record::Aux::ArrayI8(values)) if values.len() > 1 => {
                return Ok(());
            }
            Ok(_) => invalid = true,
            Err(_) => {}
        }
    }
    if checked == 0 {
        bail!("dwell time requires BAM alignments, but the input is empty");
    }
    bail!(
        "no valid mv move-table tag was found in the first {checked} alignments{}",
        if invalid {
            "; an mv tag had no values"
        } else {
            ""
        }
    );
}

pub fn run_pipeline(options: &RunOptions) -> Result<()> {
    if options.threads == 0 {
        bail!("threads must be positive");
    }
    if options.chunk_size <= 0 {
        bail!("chunk size must be positive");
    }
    if options.bed.is_some() && options.known_vcf.is_some() {
        bail!("provide either --bed_fn or --vcf_fn, not both");
    }
    if options.haploid_precise && options.haploid_sensitive {
        bail!("haploid precise and sensitive modes are mutually exclusive");
    }
    if options.enable_dwell_time && options.platform != Platform::Ont {
        bail!("dwell time is supported only for ONT");
    }
    validate_tools(options)?;
    if options.enable_dwell_time {
        validate_dwell_tags(&options.bam, &options.reference)?;
    }
    if options.use_longphase_for_intermediate_phasing && options.platform == Platform::Ilmn {
        bail!("LongPhase intermediate phasing is not supported for Illumina");
    }
    if options.gvcf && options.gq_bin_size <= 0 {
        bail!("GVCF GQ bin size must be positive");
    }
    let pileup_checkpoint = checkpoint_path(&options.model_path, &options.pileup_model_prefix)?;
    let full_checkpoint =
        checkpoint_path(&options.model_path, &options.full_alignment_model_prefix)?;
    let device = super::call::device_from_options(options.use_gpu, &options.device)?;
    let pileup_device = if matches!(device, Device::Mps | Device::Cuda(_)) {
        Device::Cpu
    } else {
        device
    };
    let reference = Reference::open(&options.reference)?;
    let contigs = prepare_contigs(&reference, &options.bam, options)?;
    fs::create_dir_all(&options.output)?;
    let tmp = options.output.join("tmp");
    let split_beds = tmp.join("split_beds");
    let phase_vcf = tmp.join("phase_output/phase_vcf");
    let phase_bam = tmp.join("phase_output/phase_bam");
    let candidate_bed = tmp.join("full_alignment_output/candidate_bed");
    let pileup_output = tmp.join("pileup_output");
    let full_output = tmp.join("full_alignment_output");
    let merge_output = tmp.join("merge_output");
    for directory in [
        &tmp,
        &split_beds,
        &phase_vcf,
        &phase_bam,
        &candidate_bed,
        &pileup_output,
        &full_output,
        &merge_output,
    ] {
        fs::create_dir_all(directory)?;
    }
    if let Some(command) = &options.command {
        fs::write(tmp.join("CMD"), command)?;
    }
    fs::write(tmp.join("CONTIGS"), contigs.join("\n"))?;
    if contigs.is_empty() {
        write_output_vcf(
            &options.output.join("merge_output.vcf"),
            &reference,
            &[],
            &contigs,
            options,
            false,
        )?;
        return Ok(());
    }
    write_extended_regions(
        &split_beds,
        options.bed.as_ref(),
        options.known_vcf.as_ref(),
        &contigs,
    )?;
    let var_pct_full = options
        .var_pct_full
        .filter(|value| *value != 0.0)
        .unwrap_or(options.platform.default_var_pct_full());
    let ref_pct_full = options
        .ref_pct_full
        .filter(|value| *value != 0.0)
        .unwrap_or(options.platform.default_ref_pct_full());
    let phasing_proportion = options
        .var_pct_phasing
        .filter(|value| *value != 0.0)
        .unwrap_or_else(|| {
            let model_name = options
                .model_path
                .file_name()
                .and_then(OsStr::to_str)
                .unwrap_or("");
            if matches!(
                model_name,
                "r941_prom_sup_g5014" | "r941_prom_hac_g5014" | "ont_guppy5"
            ) {
                0.8
            } else {
                0.7
            }
        });
    for (name, value) in [
        ("var_pct_full", var_pct_full),
        ("ref_pct_full", ref_pct_full),
        ("var_pct_phasing", phasing_proportion),
    ] {
        if !(0.0..=1.0).contains(&value) {
            bail!("{name} must be between 0 and 1");
        }
    }
    let mut jobs = Vec::new();
    let mut chunk_lines = Vec::new();
    for contig in &contigs {
        let count = if options.chunk_num > 0 {
            options.chunk_num as usize
        } else if options.chunk_num < 0 {
            1
        } else {
            ((reference.contig_len(contig) as i64 / options.chunk_size)
                + i64::from(reference.contig_len(contig) as i64 % options.chunk_size != 0))
            .max(1) as usize
        };
        for id in 1..=count {
            let mut config = call_config(
                options,
                &pileup_checkpoint,
                &options.bam,
                contig,
                pileup_device,
                var_pct_full,
            );
            config.add_indel_length = false;
            config.extend_bed = (options.bed.is_some() || options.known_vcf.is_some())
                .then(|| split_beds.join(contig));
            if options.chunk_num >= 0 {
                config.chunk_id = Some(id);
                config.chunk_num = Some(count);
            }
            config.output = Some(pileup_output.join(format!("pileup_{contig}_{id}.vcf")));
            chunk_lines.push(format!("{contig} {id} {count}"));
            jobs.push(config);
        }
    }
    fs::write(tmp.join("CHUNK_LIST"), chunk_lines.join("\n") + "\n")?;
    tch::set_num_threads(1);
    let reference_path = options.reference.clone();
    let outputs = std::thread::scope(|scope| -> Result<_> {
        let capture = (!options.no_phasing_for_fa
            && options.platform != Platform::Ilmn
            && !options.use_longphase_for_intermediate_phasing)
            .then(|| {
                scope.spawn(|| -> Result<()> {
                    for contig in &contigs {
                        crate::pipeline::phasing::capture_alignment_spans(
                            &options.bam,
                            &reference_path,
                            contig,
                            &phase_vcf.join(format!("{contig}.spans")),
                        )?;
                    }
                    Ok(())
                })
            });
        let outputs = execute_call_jobs(&jobs, options.threads.max(1), true)?;
        if let Some(capture) = capture {
            capture
                .join()
                .map_err(|_| anyhow::anyhow!("alignment span worker panicked"))??;
        }
        Ok(outputs)
    })?;
    let mut pileup = Vec::new();
    let mut nonvariants = Vec::new();
    for (config, output) in jobs.iter().zip(outputs) {
        write_vcf(
            config.output.as_deref().unwrap(),
            &reference,
            &output.records,
            &options.sample_name,
            options.command.as_deref(),
            options.gvcf,
        )?;
        pileup.extend(output.records);
        nonvariants.extend(output.nonvariants);
    }
    sort_records(&mut pileup, &contigs);
    sort_records(&mut nonvariants, &contigs);
    write_output_vcf(
        &options.output.join("pileup.vcf"),
        &reference,
        &pileup,
        &contigs,
        options,
        false,
    )?;
    if options.pileup_only {
        let filtered = merge_records(
            &pileup,
            &[],
            &MergeOptions {
                print_ref_calls: options.print_ref_calls,
                haploid_precise: options.haploid_precise,
                haploid_sensitive: options.haploid_sensitive,
                qual: options.qual,
            },
        );
        write_output_vcf(
            &options.output.join("pileup.vcf"),
            &reference,
            &filtered,
            &contigs,
            options,
            false,
        )?;
        if options.remove_intermediate_dir {
            fs::remove_dir_all(tmp)?;
        }
        return Ok(());
    }
    phase_bams(
        options,
        &reference,
        &contigs,
        &pileup,
        &read_vcf(&options.output.join("pileup.vcf.gz"))?.header,
        &phase_vcf,
        &phase_bam,
        phasing_proportion,
        var_pct_full,
    )?;
    let (variant_cutoff, reference_cutoff) = full_quality(
        &pileup,
        var_pct_full,
        Some(ref_pct_full),
        options.known_vcf.is_some(),
    );
    fs::create_dir_all(&candidate_bed)?;
    fs::write(
        candidate_bed.join("qual"),
        format!("{variant_cutoff} {reference_cutoff}"),
    )?;
    let mut full_jobs = Vec::new();
    let mut full_regions = BTreeMap::new();
    for contig in &contigs {
        let files = select_candidates(
            &options.output.join("pileup.vcf.gz"),
            None,
            &candidate_bed,
            contig,
            &options.platform.to_string(),
            &SelectOptions {
                var_pct_full,
                ref_pct_full: Some(ref_pct_full),
                qual_path: Some(candidate_bed.join("qual")),
                ..Default::default()
            },
            Some(&reference),
        )?;
        full_regions.insert(
            contig.clone(),
            candidate_bed.join(format!("full_aln_regions_{contig}")),
        );
        for file in files {
            let name = file
                .file_name()
                .and_then(OsStr::to_str)
                .context("invalid candidate filename")?;
            let bam = if options.platform == Platform::Ilmn {
                phase_bam.join(format!("{contig}.bam"))
            } else {
                options.bam.clone()
            };
            let mut config = call_config(
                options,
                &full_checkpoint,
                &bam,
                contig,
                device,
                var_pct_full,
            );
            config.full_regions = Some(file.clone());
            if options.platform != Platform::Ilmn && !options.no_phasing_for_fa {
                config.phased_vcf = Some(phase_vcf.join(format!("phased_{contig}.vcf.gz")));
            }
            config.output = Some(full_output.join(format!("full_alignment_{name}.vcf")));
            full_jobs.push(config);
        }
    }
    let full_outputs = execute_call_jobs(&full_jobs, options.threads.max(1), false)?;
    for (config, output) in full_jobs.iter().zip(full_outputs) {
        write_vcf(
            config.output.as_deref().unwrap(),
            &reference,
            &output.records,
            &options.sample_name,
            options.command.as_deref(),
            options.gvcf,
        )?;
    }
    let mut full: Vec<VcfRecord> = Vec::new();
    for path in fs::read_dir(&full_output)? {
        let path = path?.path();
        if path
            .file_name()
            .and_then(OsStr::to_str)
            .is_some_and(|name| name.starts_with("full_alignment_") && name.ends_with(".vcf"))
        {
            full.extend(read_vcf(&path)?.records);
        }
    }
    sort_records(&mut full, &contigs);
    write_output_vcf(
        &options.output.join("full_alignment.vcf"),
        &reference,
        &full,
        &contigs,
        options,
        false,
    )?;
    if options.gvcf {
        fs::create_dir_all(tmp.join("gvcf_tmp_output"))?;
        let header = vcf_header(
            &reference,
            options.command.as_deref(),
            &options.sample_name,
            true,
            &contigs,
        )?;
        write_gvcf(
            &tmp.join("gvcf_tmp_output/non_var.gvcf"),
            &header,
            &nonvariants,
        )?;
    }
    let merge_options = MergeOptions {
        print_ref_calls: options.print_ref_calls,
        haploid_precise: options.haploid_precise,
        haploid_sensitive: options.haploid_sensitive,
        qual: options.qual,
    };
    let mut merged = Vec::new();
    let mut merged_gvcf = Vec::new();
    for contig in &contigs {
        let pileup_contig: Vec<_> = pileup
            .iter()
            .filter(|record| record.contig == *contig)
            .cloned()
            .collect();
        let full_contig: Vec<_> = full
            .iter()
            .filter(|record| record.contig == *contig)
            .cloned()
            .collect();
        let records = if options.platform == Platform::Ilmn {
            let regions = IntervalMap::from_bed(&full_regions[contig], Some(contig))?;
            merge_illumina_records(
                &pileup_contig,
                &full_contig,
                &regions,
                contig,
                &merge_options,
            )
        } else {
            merge_records(&pileup_contig, &full_contig, &merge_options)
        };
        write_output_vcf(
            &merge_output.join(format!("merge_{contig}.vcf")),
            &reference,
            &records,
            &contigs,
            options,
            false,
        )?;
        if options.gvcf {
            let nonvariant_contig: Vec<_> = nonvariants
                .iter()
                .filter(|record| record.contig == *contig)
                .cloned()
                .collect();
            let gvcf_records = merge_nonvariants(
                &records,
                &nonvariant_contig,
                &reference,
                options.haploid_precise || options.haploid_sensitive,
            )?;
            write_output_vcf(
                &merge_output.join(format!("merge_{contig}.gvcf")),
                &reference,
                &gvcf_records,
                &contigs,
                options,
                true,
            )?;
            merged_gvcf.extend(gvcf_records);
        }
        merged.extend(records);
    }
    sort_records(&mut merged, &contigs);
    if let Some(path) = &options.known_vcf {
        add_back_missing(path, &mut merged)?;
        sort_records(&mut merged, &contigs);
    }
    write_output_vcf(
        &options.output.join("merge_output.vcf"),
        &reference,
        &merged,
        &contigs,
        options,
        false,
    )?;
    if options.gvcf {
        sort_records(&mut merged_gvcf, &contigs);
        write_output_vcf(
            &options.output.join("merge_output.gvcf"),
            &reference,
            &merged_gvcf,
            &contigs,
            options,
            true,
        )?;
    }
    if options.use_whatshap_for_final_phasing
        || options.use_longphase_for_final_phasing
        || options.use_whatshap_for_final_haplotagging
    {
        let mut phased_records = Vec::new();
        for contig in &contigs {
            let input = merge_output.join(format!("merge_{contig}.vcf.gz"));
            let phased = merge_output.join(format!("phased_merge_{contig}.vcf"));
            if !input.exists() || read_vcf(&input)?.records.is_empty() {
                write_output_vcf(&phased, &reference, &[], &contigs, options, false)?;
                continue;
            }
            if options.use_longphase_for_final_phasing {
                let prefix = merge_output.join(format!("phased_merge_{contig}"));
                execute(
                    Command::new(&options.longphase)
                        .arg("phase")
                        .arg("-s")
                        .arg(&input)
                        .arg("-b")
                        .arg(&options.bam)
                        .arg("-r")
                        .arg(options.reference.as_path())
                        .arg("-t")
                        .arg(options.threads.to_string())
                        .arg("-o")
                        .arg(&prefix)
                        .arg(if options.platform == Platform::Ont {
                            "--ont"
                        } else {
                            "--pb"
                        }),
                )?;
            } else {
                execute(
                    Command::new(&options.whatshap)
                        .arg("phase")
                        .arg("--output")
                        .arg(&phased)
                        .arg("--reference")
                        .arg(options.reference.as_path())
                        .arg("--ignore-read-groups")
                        .arg(&input)
                        .arg(&options.bam),
                )?;
            }
            phased_records.extend(read_vcf(&phased)?.records);
        }
        sort_records(&mut phased_records, &contigs);
        let phased_path = options.output.join("phased_merge_output.vcf");
        write_output_vcf(
            &phased_path,
            &reference,
            &phased_records,
            &contigs,
            options,
            false,
        )?;
        if options.use_whatshap_for_final_haplotagging {
            let phased_bam = options.output.join("phased_output.bam");
            execute(
                Command::new(&options.whatshap)
                    .arg("haplotag")
                    .arg("--output")
                    .arg(&phased_bam)
                    .arg("--reference")
                    .arg(options.reference.as_path())
                    .arg("--ignore-read-groups")
                    .arg(phased_path.with_extension("vcf.gz"))
                    .arg(&options.bam),
            )?;
            execute(
                Command::new(&options.samtools)
                    .arg("index")
                    .arg("-@")
                    .arg(options.threads.to_string())
                    .arg(&phased_bam),
            )?;
        }
    }
    if options.remove_intermediate_dir {
        fs::remove_dir_all(tmp)?;
    }
    Ok(())
}

pub fn check_environment(options: &RunOptions) -> Result<Vec<String>> {
    if options.threads == 0 {
        bail!("threads must be positive");
    }
    if options.chunk_size <= 0 {
        bail!("chunk size must be positive");
    }
    if options.bed.is_some() && options.known_vcf.is_some() {
        bail!("provide either --bed_fn or --vcf_fn, not both");
    }
    validate_tools(options)?;
    if options.enable_dwell_time {
        validate_dwell_tags(&options.bam, &options.reference)?;
    }
    let reference = Reference::open(&options.reference)?;
    let contigs = prepare_contigs(&reference, &options.bam, options)?;
    let tmp = options.output.join("tmp");
    for directory in [
        "log",
        "tmp/split_beds",
        "tmp/pileup_output",
        "tmp/merge_output",
        "tmp/phase_output",
        "tmp/gvcf_tmp_output",
        "tmp/full_alignment_output",
        "tmp/phase_output/phase_vcf",
        "tmp/phase_output/phase_bam",
        "tmp/full_alignment_output/candidate_bed",
    ] {
        fs::create_dir_all(options.output.join(directory))?;
    }
    fs::create_dir_all(&tmp)?;
    if let Some(command) = &options.command {
        fs::write(tmp.join("CMD"), command)?;
    }
    fs::write(tmp.join("CONTIGS"), contigs.join("\n"))?;
    if contigs.is_empty() {
        write_output_vcf(
            &options.output.join("merge_output.vcf"),
            &reference,
            &[],
            &contigs,
            options,
            false,
        )?;
        compress_index(&options.output.join("merge_output.vcf"), options)?;
        return Ok(contigs);
    }
    if options.bed.is_some() || options.known_vcf.is_some() {
        write_extended_regions(
            &tmp.join("split_beds"),
            options.bed.as_ref(),
            options.known_vcf.as_ref(),
            &contigs,
        )?;
    }
    let mut chunks = Vec::new();
    if options.chunk_num < 0 {
        chunks.push("None 0 0".to_string());
    } else {
        for contig in &contigs {
            let length = reference.contig_len(contig) as i64;
            let count = if options.chunk_num > 0 {
                options.chunk_num as usize
            } else {
                (length / options.chunk_size + i64::from(length % options.chunk_size != 0)).max(1)
                    as usize
            };
            for id in 1..=count {
                chunks.push(format!("{contig} {id} {count}"));
            }
        }
    }
    fs::write(tmp.join("CHUNK_LIST"), chunks.join("\n") + "\n")?;
    Ok(contigs)
}

pub fn check_parallel_log(path: &Path) -> Result<()> {
    let content = fs::read_to_string(path)?;
    let mut lines = content.lines();
    let Some(header) = lines.next() else {
        return Ok(());
    };
    let columns: Vec<_> = header.split('\t').collect();
    let exit = columns
        .iter()
        .position(|name| *name == "Exitval")
        .context("parallel log has no Exitval column")?;
    let signal = columns
        .iter()
        .position(|name| *name == "Signal")
        .context("parallel log has no Signal column")?;
    let command = columns
        .iter()
        .position(|name| *name == "Command")
        .context("parallel log has no Command column")?;
    for line in lines {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() <= exit.max(signal).max(command) {
            continue;
        }
        if fields[exit] != "0" || fields[signal] != "0" {
            bail!(
                "parallel command failed (exit {}, signal {}): {}",
                fields[exit],
                fields[signal],
                fields[command]
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_phasing_selection_uses_unphased_inputs() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let output = tempfile::tempdir().unwrap();
        let reference = Reference::open(&root.join("tests/fixtures/data/toy.fa")).unwrap();
        let mut options = RunOptions::new(
            root.join("tests/fixtures/data/toy.bam"),
            root.join("tests/fixtures/data/toy.fa"),
            PathBuf::new(),
            output.path().to_path_buf(),
            Platform::Ont,
        );
        options.whatshap = output
            .path()
            .join("must-not-run")
            .to_string_lossy()
            .into_owned();
        options.bgzip = root
            .join("target/legacy/env/bin/bgzip")
            .to_string_lossy()
            .into_owned();
        options.tabix = root
            .join("target/legacy/env/bin/tabix")
            .to_string_lossy()
            .into_owned();
        let header = vcf_header(
            &reference,
            None,
            "SAMPLE",
            false,
            &reference.contigs().unwrap(),
        )
        .unwrap();
        let header: Vec<_> = header.lines().map(str::to_string).collect();
        phase_bams(
            &options,
            &reference,
            &["chrToy".to_string()],
            &[],
            &header,
            &output.path().join("phase_vcf"),
            &output.path().join("phase_bam"),
            0.7,
            0.3,
        )
        .unwrap();
        let phased = read_vcf(&output.path().join("phase_vcf/phased_chrToy.vcf.gz")).unwrap();
        assert!(phased.records.is_empty());
        assert!(output.path().join("phase_bam/chrToy.bam").exists());
        assert!(output.path().join("phase_bam/chrToy.bam.bai").exists());
    }

    #[cfg(unix)]
    #[test]
    fn longphase_fallback_does_not_require_whatshap() {
        use std::os::unix::fs::PermissionsExt;

        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let output = tempfile::tempdir().unwrap();
        let phaser = output.path().join("longphase");
        fs::write(
            &phaser,
            "#!/bin/sh\nset -eu\nwhile [ \"$#\" -gt 0 ]; do\n  case \"$1\" in\n    -s) input=\"$2\"; shift 2 ;;\n    -o) output=\"$2\"; shift 2 ;;\n    *) shift ;;\n  esac\ndone\nif [ -n \"${output:-}\" ]; then cp \"$input\" \"${output}.vcf\"; fi\nexit 0\n",
        )
        .unwrap();
        let mut permissions = fs::metadata(&phaser).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&phaser, permissions).unwrap();
        let reference = Reference::open(&root.join("tests/fixtures/data/toy.fa")).unwrap();
        let mut options = RunOptions::new(
            root.join("tests/fixtures/data/toy.bam"),
            root.join("tests/fixtures/data/toy.fa"),
            PathBuf::new(),
            output.path().to_path_buf(),
            Platform::Ont,
        );
        options.whatshap = output
            .path()
            .join("missing-whatshap")
            .to_string_lossy()
            .into_owned();
        options.longphase = phaser.to_string_lossy().into_owned();
        options.bgzip = root
            .join("target/legacy/env/bin/bgzip")
            .to_string_lossy()
            .into_owned();
        options.tabix = root
            .join("target/legacy/env/bin/tabix")
            .to_string_lossy()
            .into_owned();
        let header = vcf_header(
            &reference,
            None,
            "SAMPLE",
            false,
            &reference.contigs().unwrap(),
        )
        .unwrap();
        let header: Vec<_> = header.lines().map(str::to_string).collect();
        let record = VcfRecord {
            contig: "chrToy".to_string(),
            position: 500,
            reference: "A".to_string(),
            alternate: "C".to_string(),
            quality: "20".to_string(),
            filter: "PASS".to_string(),
            id: ".".to_string(),
            info: "P".to_string(),
            format: Some("GT".to_string()),
            samples: vec!["0/1".to_string()],
        };
        phase_bams(
            &options,
            &reference,
            &["chrToy".to_string()],
            &[record],
            &header,
            &output.path().join("phase_vcf"),
            &output.path().join("phase_bam"),
            0.7,
            0.3,
        )
        .unwrap();
        let phased = read_vcf(&output.path().join("phase_vcf/phased_chrToy.vcf.gz")).unwrap();
        assert_eq!(phased.records.len(), 1);
        assert_eq!(phased.records[0].position, 500);
    }
}
