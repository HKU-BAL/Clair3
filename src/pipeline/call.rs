use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use tch::Device;

use crate::{
    call::{
        caller::{TensorInput, call_tensor_batch, write_vcf},
        decode::DecodeOptions,
    },
    config::{MAX_LONG_VARIANT_LENGTH, MAX_VARIANT_LENGTH, PREDICT_BATCH_SIZE, Platform},
    io::{
        compression::read_text_lines,
        interval::{IntervalMap, Region},
        reference::Reference,
        variant::{VcfRecord, candidate_positions},
    },
    model::{Model, ModelKind, parse_device},
    realign::workflow::{RealignmentOptions, realign_bam},
    tensor::{
        full_alignment::{
            FullAlignmentBatch, FullAlignmentOptions, PhasedVariant, candidates_from_regions,
            create_full_alignment,
        },
        pileup::{PileupBatch, PileupOptions, create_pileup},
    },
};

use super::gvcf::{GvcfOptions, nonvariant_records};

pub enum InferenceJob {
    Pileup(std::sync::Arc<PileupBatch>, DecodeOptions),
    FullAlignment(std::sync::Arc<FullAlignmentBatch>, DecodeOptions),
}

pub type InferenceRequest = (
    InferenceJob,
    std::sync::mpsc::Sender<Result<Vec<VcfRecord>>>,
);

pub enum Inference<'a> {
    Local(&'a std::sync::Mutex<Model>),
    Remote(&'a std::sync::mpsc::Sender<InferenceRequest>),
}

pub(crate) fn inference_job(job: InferenceJob, model: &Model) -> Result<Vec<VcfRecord>> {
    tch::no_grad(|| match job {
        InferenceJob::Pileup(batch, options) => call_tensor_batch(
            model,
            TensorInput {
                shape: &batch.shape,
                values: &batch.values,
                positions: &batch.positions,
                alt_info: &batch.alt_info,
            },
            &options,
            PREDICT_BATCH_SIZE,
        ),
        InferenceJob::FullAlignment(batch, options) => call_tensor_batch(
            model,
            TensorInput {
                shape: &batch.shape,
                values: &batch.values,
                positions: &batch.positions,
                alt_info: &batch.alt_info,
            },
            &options,
            PREDICT_BATCH_SIZE,
        ),
    })
}

fn run_inference(inference: &Inference, job: InferenceJob) -> Result<Vec<VcfRecord>> {
    match inference {
        Inference::Local(model) => {
            let model = model
                .lock()
                .map_err(|_| anyhow::anyhow!("inference model lock poisoned"))?;
            inference_job(job, &model)
        }
        Inference::Remote(sender) => {
            let (response, result) = std::sync::mpsc::channel();
            sender.send((job, response))?;
            result.recv().context("inference worker stopped")?
        }
    }
}

#[derive(Debug, Clone)]
pub struct CallConfig {
    pub bam: PathBuf,
    pub reference: PathBuf,
    pub checkpoint: PathBuf,
    pub output: Option<PathBuf>,
    pub platform: Platform,
    pub sample_name: String,
    pub contig: String,
    pub start: Option<i64>,
    pub end: Option<i64>,
    pub bed: Option<PathBuf>,
    pub extend_bed: Option<PathBuf>,
    pub known_vcf: Option<PathBuf>,
    pub full_regions: Option<PathBuf>,
    pub min_coverage: i64,
    pub min_mq: u8,
    pub min_bq: u8,
    pub min_af: Option<f64>,
    pub snp_min_af: f64,
    pub indel_min_af: f64,
    pub fast_mode: bool,
    pub call_snp_only: bool,
    pub call_head_tail: bool,
    pub enable_long_indel: bool,
    pub enable_dwell_time: bool,
    pub add_indel_length: bool,
    pub haploid_precise: bool,
    pub haploid_sensitive: bool,
    pub no_phasing_for_fa: bool,
    pub add_no_phasing_data_training: bool,
    pub phasing_info_in_bam: bool,
    pub need_realignment: bool,
    pub gvcf: bool,
    pub gvcf_options: GvcfOptions,
    pub chunk_id: Option<usize>,
    pub chunk_num: Option<usize>,
    pub device: Device,
    pub show_reference: bool,
    pub quality_for_pass: Option<f64>,
    pub keep_iupac: bool,
    pub phased_vcf: Option<PathBuf>,
    pub command: Option<String>,
}

impl CallConfig {
    pub fn new(
        bam: PathBuf,
        reference: PathBuf,
        checkpoint: PathBuf,
        contig: String,
        platform: Platform,
    ) -> Self {
        Self {
            bam,
            reference,
            checkpoint,
            output: None,
            platform,
            sample_name: "SAMPLE".to_string(),
            contig,
            start: None,
            end: None,
            bed: None,
            extend_bed: None,
            known_vcf: None,
            full_regions: None,
            min_coverage: 2,
            min_mq: 5,
            min_bq: 0,
            min_af: None,
            snp_min_af: 0.08,
            indel_min_af: platform.min_indel_af(),
            fast_mode: false,
            call_snp_only: false,
            call_head_tail: false,
            enable_long_indel: false,
            enable_dwell_time: false,
            add_indel_length: true,
            haploid_precise: false,
            haploid_sensitive: false,
            no_phasing_for_fa: false,
            add_no_phasing_data_training: false,
            phasing_info_in_bam: false,
            need_realignment: true,
            gvcf: false,
            gvcf_options: GvcfOptions::default(),
            chunk_id: None,
            chunk_num: None,
            device: Device::Cpu,
            show_reference: true,
            quality_for_pass: Some(2.0),
            keep_iupac: false,
            phased_vcf: None,
            command: None,
        }
    }
}

#[derive(Debug, Default)]
pub struct CallOutput {
    pub records: Vec<VcfRecord>,
    pub nonvariants: Vec<VcfRecord>,
}

pub fn calling_region(reference: &Reference, config: &CallConfig) -> Result<Region> {
    let contig_length = reference.contig_len(&config.contig) as i64;
    if contig_length == 0 {
        bail!("contig {} is absent from the reference", config.contig);
    }
    if config.chunk_id.is_some() != config.chunk_num.is_some() {
        bail!("chunk ID and chunk count must be supplied together");
    }
    if let (Some(id), Some(count)) = (config.chunk_id, config.chunk_num)
        && (id == 0 || count == 0 || id > count)
    {
        bail!("chunk ID must be between 1 and the chunk count");
    }
    if let Some(path) = &config.known_vcf {
        let mut positions = candidate_positions(path, Some(&config.contig))?;
        if let (Some(id), Some(count)) = (config.chunk_id, config.chunk_num) {
            let chunk_size = positions.len().div_ceil(count);
            let start = (id - 1) * chunk_size;
            positions = positions.into_iter().skip(start).take(chunk_size).collect();
        }
        let (Some(first), Some(last)) = (positions.first(), positions.last()) else {
            return Region::new(&config.contig, 0, 0);
        };
        return Region::new(&config.contig, first - 1, *last);
    }
    let (start, end) = if let Some(path) = &config.extend_bed {
        let regions = IntervalMap::from_bed(path, Some(&config.contig))?;
        let (bed_start, bed_end) = regions.bounds(&config.contig).with_context(|| {
            format!("contig {} is absent from {}", config.contig, path.display())
        })?;
        (bed_start, bed_end)
    } else if let (Some(start), Some(end)) = (config.start, config.end) {
        if start < 1 || end < start || end > contig_length {
            bail!(
                "invalid one-based contig range {}:{}-{end}",
                config.contig,
                start
            );
        }
        (start - 1, end)
    } else {
        (0, contig_length)
    };
    if let (Some(id), Some(count)) = (config.chunk_id, config.chunk_num) {
        let length = end - start;
        let size = length / count as i64 + i64::from(length % count as i64 != 0);
        let start = start + size * (id - 1) as i64;
        return Region::new(&config.contig, start, (start + size).min(end));
    }
    Region::new(&config.contig, start, end)
}

fn known_positions(config: &CallConfig) -> Result<Option<BTreeSet<i64>>> {
    config
        .known_vcf
        .as_ref()
        .map(|path| {
            candidate_positions(path, Some(&config.contig))
                .map(|positions| positions.into_iter().collect())
        })
        .transpose()
}

pub(crate) fn pileup_batch(
    config: &CallConfig,
    reference: &Reference,
) -> Result<(PileupBatch, Region)> {
    let region = calling_region(reference, config)?;
    if region.start == region.end {
        return Ok((PileupBatch::default(), region));
    }
    let confident_regions = config
        .extend_bed
        .as_ref()
        .or(config.bed.as_ref())
        .map(|path| IntervalMap::from_bed(path, Some(&config.contig)))
        .transpose()?;
    let known_positions = known_positions(config)?;
    let expanded = Region::new(
        &config.contig,
        (region.start - 33).max(0),
        (region.end + 33).min(reference.contig_len(&config.contig) as i64),
    )?;
    let options = PileupOptions {
        min_depth: config.min_coverage,
        illumina: config.platform == Platform::Ilmn,
        min_snp_af: config.min_af.unwrap_or(config.snp_min_af),
        min_indel_af: config.min_af.unwrap_or(config.indel_min_af),
        min_mq: config.min_mq,
        max_indel_length: if config.enable_long_indel {
            MAX_LONG_VARIANT_LENGTH
        } else {
            MAX_VARIANT_LENGTH
        },
        call_snp_only: config.call_snp_only,
        call_head_tail: config.call_head_tail,
        fast_mode: config.fast_mode && config.platform == Platform::Ont,
        gvcf: config.gvcf,
        confident_regions: confident_regions.as_ref(),
        known_positions: known_positions.as_ref(),
        ..Default::default()
    };
    Ok((
        create_pileup(&config.bam, reference, &expanded, &options)?,
        region,
    ))
}

pub(crate) fn call_pileup(
    config: &CallConfig,
    reference: &Reference,
    inference: &Inference,
) -> Result<CallOutput> {
    let profile_start = std::time::Instant::now();
    let (batch, region) = pileup_batch(config, reference)?;
    let tensor_elapsed = profile_start.elapsed();
    let inference_start = std::time::Instant::now();
    if region.start == region.end {
        return Ok(Default::default());
    }
    let batch = std::sync::Arc::new(batch);
    let records = if batch.positions.is_empty() {
        Vec::new()
    } else {
        run_inference(
            inference,
            InferenceJob::Pileup(
                batch.clone(),
                DecodeOptions {
                    show_reference: config.show_reference,
                    haploid_precise: config.haploid_precise,
                    haploid_sensitive: config.haploid_sensitive,
                    quality_for_pass: config.quality_for_pass,
                    add_indel_length: config.add_indel_length,
                    pileup: true,
                    gvcf: config.gvcf,
                    enable_long_indel: config.enable_long_indel,
                    maximum_variant_length: if config.enable_long_indel {
                        MAX_LONG_VARIANT_LENGTH
                    } else {
                        MAX_VARIANT_LENGTH
                    },
                    keep_iupac: config.keep_iupac,
                },
            ),
        )?
    };
    if std::env::var_os("CLAIR3_PROFILE").is_some() {
        eprintln!(
            "pileup {} tensor={:.3}s inference={:.3}s",
            config.contig,
            tensor_elapsed.as_secs_f64(),
            inference_start.elapsed().as_secs_f64()
        );
    }
    let records = records
        .into_iter()
        .filter(|record| {
            record.contig == config.contig
                && record.position > region.start
                && record.position <= region.end
        })
        .collect();
    let nonvariants = if config.gvcf {
        let expanded = Region::new(
            &region.contig,
            (region.start - 33).max(0),
            (region.end + 33).min(reference.contig_len(&region.contig) as i64),
        )?;
        nonvariant_records(reference, &expanded, &batch, &config.gvcf_options)?
            .into_iter()
            .filter(|record| record.position > region.start && record.position < region.end)
            .collect()
    } else {
        Vec::new()
    };
    Ok(CallOutput {
        records,
        nonvariants,
    })
}

fn phased_from_regions(path: &Path, contig: &str) -> Result<Vec<PhasedVariant>> {
    let mut output = Vec::new();
    for line in read_text_lines(path)? {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 4 || fields[0] != contig {
            continue;
        }
        let metadata: Vec<_> = fields[3].split('-').collect();
        if metadata.len() != 4 {
            bail!("invalid phased variant metadata in {}", path.display());
        }
        output.push(PhasedVariant {
            position: fields[1].parse::<i64>()? + 11,
            ref_base: metadata[0].as_bytes()[0],
            alt_base: metadata[1].as_bytes()[0],
            genotype: metadata[2].parse()?,
            phase_set: metadata[3].parse()?,
        });
    }
    output.sort_by_key(|variant| variant.position);
    Ok(output)
}

fn full_region(
    path: &Path,
    contig: &str,
    reference: &Reference,
    platform: Platform,
) -> Result<Option<Region>> {
    let mut start = i64::MAX;
    let mut end = 0_i64;
    for line in read_text_lines(path)? {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 3 || fields[0] != contig {
            continue;
        }
        start = start.min(fields[1].parse::<i64>()?);
        end = end.max(fields[2].parse::<i64>()?);
    }
    if start == i64::MAX {
        return Ok(None);
    }
    let length = reference.contig_len(contig) as i64;
    let (start, end) = if platform == Platform::Ilmn {
        ((start + 1 - 1000).max(0), (end + 1001).min(length))
    } else {
        (start.max(0), end.min(length))
    };
    Ok(Some(Region::new(contig, start, end)?))
}

pub(crate) fn full_alignment_batch(
    config: &CallConfig,
    reference: &Reference,
) -> Result<FullAlignmentBatch> {
    let (region, mut candidates, mut variants) = if let Some(path) = &config.full_regions {
        let Some(region) = full_region(path, &config.contig, reference, config.platform)? else {
            return Ok(Default::default());
        };
        (
            region,
            candidates_from_regions(path, &config.contig)?,
            phased_from_regions(path, &config.contig)?,
        )
    } else {
        (calling_region(reference, config)?, Vec::new(), Vec::new())
    };
    if variants.is_empty()
        && !config.no_phasing_for_fa
        && let Some(path) = &config.phased_vcf
    {
        variants = crate::tensor::full_alignment::phased_variants_from_vcf(path, &config.contig)?;
    }
    if region.start == region.end {
        return Ok(Default::default());
    }
    let realigned;
    let bam = if config.platform == Platform::Ilmn && config.need_realignment {
        let name = config
            .bam
            .file_name()
            .context("input BAM has no filename")?
            .to_string_lossy();
        realigned = config
            .bam
            .with_file_name(format!("{name}.{}_{}", region.start, region.end));
        realign_bam(
            &config.bam,
            reference,
            &config.contig,
            Some(region.start + 1),
            Some(region.end),
            &realigned,
            &RealignmentOptions {
                candidate_bed: config.full_regions.clone(),
                extend_bed: config.extend_bed.clone(),
                min_coverage: config.min_coverage as f64,
                min_mq: config.min_mq,
                chunk_id: config.chunk_id,
                chunk_num: config.chunk_num,
            },
        )?;
        realigned.as_path()
    } else {
        config.bam.as_path()
    };
    if config.platform == Platform::Ilmn || candidates.is_empty() {
        let candidate_pileup = create_pileup(
            bam,
            reference,
            &region,
            &PileupOptions {
                min_depth: config.min_coverage,
                illumina: config.platform == Platform::Ilmn,
                min_snp_af: config.min_af.unwrap_or(config.snp_min_af),
                min_indel_af: if config.platform == Platform::Ilmn {
                    config.min_af.unwrap_or(0.15)
                } else {
                    config.min_af.unwrap_or(config.indel_min_af)
                },
                min_mq: config.min_mq,
                max_indel_length: if config.enable_long_indel {
                    MAX_LONG_VARIANT_LENGTH
                } else {
                    MAX_VARIANT_LENGTH
                },
                exclude_deletions_from_af_denominator: config.platform == Platform::Ilmn,
                allow_empty_flanking: config.platform == Platform::Ilmn,
                call_head_tail: config.platform == Platform::Ilmn || config.call_head_tail,
                ..Default::default()
            },
        )?;
        candidates = candidate_pileup
            .positions
            .iter()
            .map(|position| {
                position
                    .rsplit(':')
                    .nth(1)
                    .context("invalid candidate position")?
                    .parse::<i64>()
                    .map(|position| position - 1)
                    .map_err(anyhow::Error::from)
            })
            .collect::<Result<_>>()?;
    }
    if let Some(positions) = known_positions(config)? {
        candidates.retain(|candidate| positions.contains(&(candidate + 1)));
    }
    candidates.sort_unstable();
    candidates.dedup();
    if candidates.is_empty() {
        return Ok(Default::default());
    }
    let options = FullAlignmentOptions {
        matrix_depth: config.platform.matrix_depth(),
        max_indel_length: if config.enable_long_indel {
            MAX_LONG_VARIANT_LENGTH
        } else {
            MAX_VARIANT_LENGTH
        },
        min_mq: config.min_mq,
        min_bq: config.min_bq,
        need_haplotagging: !config.no_phasing_for_fa,
        phasing_info_in_bam: config.phasing_info_in_bam,
        enable_dwell_time: config.enable_dwell_time,
        illumina: config.platform == Platform::Ilmn,
        ..Default::default()
    };
    let mut batch =
        create_full_alignment(bam, reference, &region, &candidates, &variants, &options)?;
    if config.add_no_phasing_data_training {
        let augmented = create_full_alignment(
            bam,
            reference,
            &region,
            &candidates,
            &variants,
            &FullAlignmentOptions {
                need_haplotagging: false,
                phasing_info_in_bam: false,
                force_unphased: true,
                ..options
            },
        )?;
        batch.append(augmented)?;
    }
    Ok(batch)
}

pub(crate) fn call_full_alignment(
    config: &CallConfig,
    reference: &Reference,
    inference: &Inference,
) -> Result<CallOutput> {
    let profile_start = std::time::Instant::now();
    let batch = full_alignment_batch(config, reference)?;
    let tensor_elapsed = profile_start.elapsed();
    let inference_start = std::time::Instant::now();
    if batch.positions.is_empty() {
        return Ok(Default::default());
    }
    let records = run_inference(
        inference,
        InferenceJob::FullAlignment(
            std::sync::Arc::new(batch),
            DecodeOptions {
                show_reference: config.show_reference,
                haploid_precise: config.haploid_precise,
                haploid_sensitive: config.haploid_sensitive,
                quality_for_pass: config.quality_for_pass,
                add_indel_length: config.add_indel_length,
                gvcf: config.gvcf,
                enable_long_indel: config.enable_long_indel,
                maximum_variant_length: if config.enable_long_indel {
                    MAX_LONG_VARIANT_LENGTH
                } else {
                    MAX_VARIANT_LENGTH
                },
                keep_iupac: config.keep_iupac,
                ..Default::default()
            },
        ),
    )?;
    if std::env::var_os("CLAIR3_PROFILE").is_some() {
        eprintln!(
            "full {} tensor={:.3}s inference={:.3}s",
            config.contig,
            tensor_elapsed.as_secs_f64(),
            inference_start.elapsed().as_secs_f64()
        );
    }
    Ok(CallOutput {
        records,
        nonvariants: Vec::new(),
    })
}

pub fn call_var_bam(config: &CallConfig, pileup: bool) -> Result<CallOutput> {
    if config.checkpoint.as_os_str().is_empty() {
        bail!("a model checkpoint is required for variant calling");
    }
    if config.enable_dwell_time && (pileup || config.platform != Platform::Ont) {
        bail!("dwell time is available only for ONT full-alignment calling");
    }
    let reference = Reference::open(&config.reference)?;
    let model = Model::new(
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
        config.add_indel_length,
        Some(&config.checkpoint),
    )?;
    let model = std::sync::Mutex::new(model);
    let inference = Inference::Local(&model);
    let output = if pileup {
        call_pileup(config, &reference, &inference)?
    } else {
        call_full_alignment(config, &reference, &inference)?
    };
    if let Some(path) = &config.output {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        write_vcf(
            path,
            &reference,
            &output.records,
            &config.sample_name,
            config.command.as_deref(),
            config.gvcf,
        )?;
    }
    Ok(output)
}

pub fn checkpoint_path(model_path: &Path, prefix: &str) -> Result<PathBuf> {
    let path = model_path.join(prefix);
    if path.extension().and_then(|extension| extension.to_str()) == Some("pt") {
        return Ok(path);
    }
    let with_extension = path.with_extension("pt");
    if with_extension.is_file() {
        return Ok(with_extension);
    }
    if path.is_file() {
        return Ok(path);
    }
    bail!("checkpoint {} does not exist", with_extension.display())
}

pub fn device_from_options(use_gpu: bool, device: &str) -> Result<Device> {
    if device != "EMPTY" && !device.is_empty() {
        return parse_device(device);
    }
    if use_gpu {
        if !tch::Cuda::is_available() {
            bail!("GPU calling was requested, but CUDA is unavailable");
        }
        Ok(Device::Cuda(0))
    } else {
        Ok(Device::Cpu)
    }
}
