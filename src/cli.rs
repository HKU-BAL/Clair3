use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::{ArgAction, Args, CommandFactory, Parser, Subcommand};

use crate::{
    VERSION,
    config::Platform,
    io::{
        reference::Reference,
        variant::{VcfRecord, vcf_header},
    },
    model::{Model, ModelKind, parse_device},
    pipeline::{
        call::{CallConfig, call_var_bam, checkpoint_path, device_from_options},
        gvcf::{GvcfOptions, merge_nonvariants, write_gvcf},
        merge::{MergeOptions, merge_illumina_records, merge_records, read_vcf},
        run::{RunOptions, check_parallel_log, run_pipeline},
        select::{SelectOptions, full_quality, phase_quality, select_candidates, select_het_snp},
    },
    realign::workflow::{RealignmentOptions, realign_bam},
};

#[derive(Parser)]
#[command(name = "clair3", version = VERSION, about = "Clair3 variant calling and model training in Rust")]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    #[command(name = "run", about = "Call variants from an indexed BAM/CRAM")]
    Run(RunArgs),
    #[command(name = "CallVarBam")]
    CallVarBam(CallVarBamArgs),
    #[command(name = "RealignReads")]
    RealignReads(RealignArgs),
    #[command(name = "CallVariants")]
    CallVariants(CallVariantsArgs),
    #[command(name = "CreateTensorPileup")]
    CreateTensorPileup(CallVarBamArgs),
    #[command(name = "CreateTensorFullAlignment")]
    CreateTensorFullAlignment(CallVarBamArgs),
    #[command(name = "CreateTrainingTensor")]
    CreateTrainingTensor(TrainingTensorArgs),
    #[command(name = "GetTruth")]
    GetTruth(crate::pipeline::training::GetTruthArgs),
    #[command(name = "SplitExtendBed")]
    SplitExtendBed(crate::pipeline::training::SplitExtendBedArgs),
    #[command(name = "UnifyRepresentation")]
    UnifyRepresentation(crate::pipeline::unification::UnifyRepresentationArgs),
    #[command(name = "SelectQual")]
    SelectQual(SelectQualArgs),
    #[command(name = "SelectHetSnp")]
    SelectHetSnp(SelectHetSnpArgs),
    #[command(name = "SelectCandidates")]
    SelectCandidates(SelectCandidatesArgs),
    #[command(name = "MergeVcf")]
    MergeVcf(MergeVcfArgs),
    #[command(name = "CheckExitCode")]
    CheckExitCode(CheckExitCodeArgs),
    #[command(name = "CheckEnvs")]
    CheckEnvs(crate::pipeline::commands::CheckEnvsArgs),
    #[command(name = "SortVcf")]
    SortVcf(crate::pipeline::commands::SortVcfArgs),
    #[command(name = "AddBackMissingVariantsInGenotyping")]
    AddBackMissingVariantsInGenotyping(crate::pipeline::commands::AddBackMissingArgs),
    #[command(name = "AddPairEndAlleleDepth")]
    AddPairEndAlleleDepth(crate::pipeline::commands::PairEndDepthArgs),
    #[command(name = "SwitchZygosityBasedOnSVCalls")]
    SwitchZygosityBasedOnSVCalls(crate::pipeline::commands::SwitchZygosityArgs),
    #[command(name = "GetOverallMetrics")]
    GetOverallMetrics(crate::pipeline::commands::OverallMetricsArgs),
    #[command(name = "Train", visible_alias = "train")]
    Train(TrainArgs),
    #[command(name = "Tensor2Bin")]
    Tensor2Bin(TensorBinArgs),
    #[command(name = "MergeBin")]
    MergeBin(MergeBinArgs),
    #[command(
        name = "ConvertCheckpoint",
        visible_alias = "ConvertTfCheckpoint",
        visible_alias = "convert-tf-checkpoint"
    )]
    ConvertCheckpoint(ConvertCheckpointArgs),
    #[command(name = "model-info", visible_alias = "ModelInfo")]
    ModelInfo(ModelInfoArgs),
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct CommonCallArgs {
    #[arg(short = 'p', long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long = "min_mq", alias = "minMQ", default_value_t = 5)]
    min_mq: u8,
    #[arg(long = "min_bq", alias = "minBQ", default_value_t = 0)]
    min_bq: u8,
    #[arg(long = "min_coverage", alias = "minCoverage", default_value_t = 2)]
    min_coverage: i64,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    fast_mode: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    call_snp_only: bool,
    #[arg(long = "enable_variant_calling_at_sequence_head_and_tail", num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    call_head_tail: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    enable_long_indel: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    keep_iupac_bases: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    enable_dwell_time: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    no_phasing_for_fa: bool,
    #[arg(long = "need_realignment", action = ArgAction::SetFalse, default_value_t = true, hide = true)]
    need_realignment: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    gvcf: bool,
    #[arg(long, default_value_t = 0.001)]
    base_err: f64,
    #[arg(long, default_value_t = 5)]
    gq_bin_size: i32,
    #[arg(long, default_value_t = false)]
    bp_resolution: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    print_ref_calls: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_precise: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_sensitive: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_gpu: bool,
    #[arg(long, default_value = "EMPTY")]
    device: String,
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct RunArgs {
    #[arg(short = 'b', long = "bam_fn")]
    bam: PathBuf,
    #[arg(short = 'f', long = "ref_fn")]
    reference: PathBuf,
    #[arg(short = 'm', long = "model_path")]
    model_path: PathBuf,
    #[arg(short = 'o', long = "output")]
    output: PathBuf,
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,
    #[arg(long = "bed_fn", default_value = "EMPTY")]
    bed: String,
    #[arg(long = "vcf_fn", default_value = "EMPTY")]
    known_vcf: String,
    #[arg(long = "ctg_name", default_value = "EMPTY")]
    contig_names: String,
    #[arg(long = "sample_name", alias = "sampleName", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long, default_value_t = 0)]
    chunk_num: isize,
    #[arg(long, default_value_t = 5_000_000)]
    chunk_size: i64,
    #[arg(long, default_value_t = 2.0)]
    qual: f64,
    #[arg(long, default_value = "samtools")]
    samtools: String,
    #[arg(long, default_value = "whatshap")]
    whatshap: String,
    #[arg(long, default_value = "EMPTY")]
    longphase: String,
    #[arg(long, default_value = "bgzip", hide = true)]
    bgzip: String,
    #[arg(long, default_value = "tabix", hide = true)]
    tabix: String,
    #[arg(long, default_value_t = 0.0)]
    var_pct_full: f64,
    #[arg(long, default_value_t = 0.0)]
    ref_pct_full: f64,
    #[arg(long, default_value_t = 0.0)]
    var_pct_phasing: f64,
    #[arg(long, default_value_t = 0.0)]
    snp_min_af: f64,
    #[arg(long, default_value_t = 0.0)]
    indel_min_af: f64,
    #[arg(long, default_value_t = 0)]
    min_contig_size: u64,
    #[arg(long, default_value = "pileup")]
    pileup_model_prefix: String,
    #[arg(long = "fa_model_prefix", default_value = "full_alignment")]
    full_alignment_model_prefix: String,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    pileup_only: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    output_all_contigs_in_gvcf_header: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    include_all_ctgs: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    remove_intermediate_dir: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    enable_phasing: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    longphase_for_phasing: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_longphase_for_intermediate_phasing: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = true, hide = true)]
    use_whatshap_for_intermediate_phasing: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_whatshap_for_final_output_phasing: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_longphase_for_final_output_phasing: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_whatshap_for_final_output_haplotagging: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false, hide = true)]
    disable_c_impl: bool,
    #[arg(long, hide = true)]
    python: Option<String>,
    #[arg(long, hide = true)]
    pypy: Option<String>,
    #[arg(long, hide = true)]
    parallel: Option<String>,
    #[command(flatten)]
    common: CommonCallArgs,
}

fn empty_path(value: &str) -> Option<PathBuf> {
    (value != "EMPTY" && !value.is_empty()).then(|| PathBuf::from(value))
}

impl RunArgs {
    fn options(&self) -> Result<RunOptions> {
        let _ = (
            self.use_whatshap_for_intermediate_phasing,
            self.disable_c_impl,
            &self.python,
            &self.pypy,
            &self.parallel,
        );
        let mut options = RunOptions::new(
            self.bam.clone(),
            self.reference.clone(),
            self.model_path.clone(),
            self.output.clone(),
            self.common.platform,
        );
        options.threads = self.threads;
        options.bed = empty_path(&self.bed);
        options.known_vcf = empty_path(&self.known_vcf);
        options.contig_names = (self.contig_names != "EMPTY").then(|| self.contig_names.clone());
        options.sample_name = self.sample_name.clone();
        options.chunk_num = self.chunk_num;
        options.chunk_size = self.chunk_size;
        options.qual = (self.qual >= 0.0).then_some(self.qual);
        options.samtools = self.samtools.clone();
        options.whatshap = self.whatshap.clone();
        options.longphase = if self.longphase == "EMPTY" {
            "longphase".to_string()
        } else {
            self.longphase.clone()
        };
        options.bgzip = self.bgzip.clone();
        options.tabix = self.tabix.clone();
        options.var_pct_full = (self.var_pct_full != 0.0).then_some(self.var_pct_full);
        options.ref_pct_full = (self.ref_pct_full != 0.0).then_some(self.ref_pct_full);
        options.var_pct_phasing = (self.var_pct_phasing != 0.0).then_some(self.var_pct_phasing);
        options.snp_min_af =
            (self.snp_min_af != 0.0 || options.known_vcf.is_some()).then_some(self.snp_min_af);
        options.indel_min_af =
            (self.indel_min_af != 0.0 || options.known_vcf.is_some()).then_some(self.indel_min_af);
        options.min_mq = self.common.min_mq.max(5);
        options.min_bq = self.common.min_bq;
        options.min_coverage = self.common.min_coverage.max(2);
        options.min_contig_size = self.min_contig_size;
        options.pileup_model_prefix = self.pileup_model_prefix.clone();
        options.full_alignment_model_prefix = self.full_alignment_model_prefix.clone();
        options.pileup_only = self.pileup_only;
        options.output_all_contigs_in_gvcf_header = self.output_all_contigs_in_gvcf_header;
        options.include_all_contigs = self.include_all_ctgs;
        options.remove_intermediate_dir = self.remove_intermediate_dir;
        options.use_longphase_for_intermediate_phasing =
            self.use_longphase_for_intermediate_phasing || self.longphase_for_phasing;
        options.use_whatshap_for_final_phasing = self.use_whatshap_for_final_output_phasing
            || self.enable_phasing
            || self.use_whatshap_for_final_output_haplotagging;
        options.use_longphase_for_final_phasing = self.use_longphase_for_final_output_phasing;
        options.use_whatshap_for_final_haplotagging =
            self.use_whatshap_for_final_output_haplotagging;
        options.fast_mode = self.common.fast_mode;
        options.call_snp_only = self.common.call_snp_only;
        options.call_head_tail = self.common.call_head_tail;
        options.enable_long_indel = self.common.enable_long_indel;
        options.keep_iupac = self.common.keep_iupac_bases;
        options.enable_dwell_time = self.common.enable_dwell_time;
        options.no_phasing_for_fa = self.common.no_phasing_for_fa;
        options.need_realignment = self.common.need_realignment;
        options.gvcf = self.common.gvcf;
        options.gq_bin_size = self.common.gq_bin_size;
        options.base_error = self.common.base_err;
        options.bp_resolution = self.common.bp_resolution;
        options.print_ref_calls = self.common.print_ref_calls;
        options.haploid_precise = self.common.haploid_precise;
        options.haploid_sensitive = self.common.haploid_sensitive;
        options.use_gpu = self.common.use_gpu;
        options.device = self.common.device.clone();
        options.command = Some(std::env::args().collect::<Vec<_>>().join(" "));
        Ok(options)
    }
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct CallVarBamArgs {
    #[arg(long = "bam_fn")]
    bam: PathBuf,
    #[arg(long = "chkpnt_fn")]
    checkpoint: Option<PathBuf>,
    #[arg(long = "tensor_can_fn")]
    tensor_output: Option<PathBuf>,
    #[arg(long = "ref_fn")]
    reference: PathBuf,
    #[arg(long = "call_fn")]
    output: Option<PathBuf>,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long = "ctgStart")]
    start: Option<i64>,
    #[arg(long = "ctgEnd")]
    end: Option<i64>,
    #[arg(long = "bed_fn")]
    bed: Option<PathBuf>,
    #[arg(long = "extend_bed")]
    extend_bed: Option<PathBuf>,
    #[arg(long = "vcf_fn")]
    known_vcf: Option<PathBuf>,
    #[arg(long = "full_aln_regions")]
    full_regions: Option<PathBuf>,
    #[arg(long = "phased_vcf_fn")]
    phased_vcf: Option<PathBuf>,
    #[arg(long = "sampleName", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long)]
    min_af: Option<f64>,
    #[arg(long, default_value_t = 0.08)]
    snp_min_af: f64,
    #[arg(long)]
    indel_min_af: Option<f64>,
    #[arg(long, default_value_t = false)]
    pileup: bool,
    #[arg(long, default_value_t = false)]
    add_indel_length: bool,
    #[arg(long = "showRef", action = ArgAction::SetFalse, default_value_t = true, hide = true)]
    show_reference: bool,
    #[arg(long = "phasing_info_in_bam", default_value_t = false, hide = true)]
    phasing_info_in_bam: bool,
    #[arg(long = "unify_repre", default_value_t = false, hide = true)]
    unify_representation: bool,
    #[arg(long = "unify_repre_fn", hide = true)]
    unify_output: Option<PathBuf>,
    #[arg(long)]
    chunk_id: Option<usize>,
    #[arg(long)]
    chunk_num: Option<usize>,
    #[arg(long = "temp_file_dir", default_value = "./")]
    temp_file_dir: PathBuf,
    #[arg(long = "cmd_fn")]
    command_file: Option<PathBuf>,
    #[arg(long = "samtools", hide = true)]
    samtools: Option<String>,
    #[arg(long = "python", hide = true)]
    python: Option<String>,
    #[arg(long = "pypy", hide = true)]
    pypy: Option<String>,
    #[arg(long = "delay", hide = true)]
    delay: Option<usize>,
    #[arg(
        long = "torch_threads",
        alias = "tensorflow_threads",
        default_value_t = 4
    )]
    torch_threads: i32,
    #[command(flatten)]
    common: CommonCallArgs,
}

fn resolve_checkpoint(path: &Path) -> Result<PathBuf> {
    if path.extension().and_then(|extension| extension.to_str()) == Some("pt") {
        return Ok(path.to_path_buf());
    }
    checkpoint_path(
        path.parent().unwrap_or_else(|| Path::new(".")),
        &path
            .file_name()
            .context("checkpoint path has no filename")?
            .to_string_lossy(),
    )
}

fn call_config(args: &CallVarBamArgs) -> Result<CallConfig> {
    if args.torch_threads <= 0 {
        bail!("torch_threads must be positive");
    }
    tch::set_num_threads(args.torch_threads);
    let checkpoint = args
        .checkpoint
        .as_ref()
        .map(|path| resolve_checkpoint(path))
        .transpose()?
        .unwrap_or_default();
    let mut config = CallConfig::new(
        args.bam.clone(),
        args.reference.clone(),
        checkpoint,
        args.contig.clone(),
        args.common.platform,
    );
    config.output = args.output.clone();
    config.sample_name = args.sample_name.clone();
    config.start = args.start;
    config.end = args.end;
    config.bed = args.bed.clone();
    config.extend_bed = args.extend_bed.clone();
    config.known_vcf = args.known_vcf.clone();
    config.full_regions = args.full_regions.clone();
    config.phased_vcf = args.phased_vcf.clone();
    config.min_mq = args.common.min_mq;
    config.min_bq = args.common.min_bq;
    config.min_coverage = args.common.min_coverage;
    config.min_af = args.min_af;
    config.snp_min_af = args.snp_min_af;
    config.indel_min_af = args
        .indel_min_af
        .unwrap_or(args.common.platform.min_indel_af());
    config.fast_mode = args.common.fast_mode;
    config.call_snp_only = args.common.call_snp_only;
    config.call_head_tail = args.common.call_head_tail;
    config.enable_long_indel = args.common.enable_long_indel;
    config.keep_iupac = args.common.keep_iupac_bases;
    config.enable_dwell_time = args.common.enable_dwell_time;
    config.add_indel_length = args.add_indel_length;
    config.haploid_precise = args.common.haploid_precise;
    config.haploid_sensitive = args.common.haploid_sensitive;
    config.no_phasing_for_fa = args.common.no_phasing_for_fa;
    config.phasing_info_in_bam = args.phasing_info_in_bam;
    config.need_realignment = args.common.need_realignment;
    config.gvcf = args.common.gvcf;
    config.gvcf_options = GvcfOptions {
        base_error: args.common.base_err,
        gq_bin_size: args.common.gq_bin_size,
        bp_resolution: args.common.bp_resolution,
    };
    config.chunk_id = args.chunk_id;
    config.chunk_num = args.chunk_num;
    config.device = device_from_options(args.common.use_gpu, &args.common.device)?;
    config.show_reference = args.show_reference;
    config.command = args
        .command_file
        .as_ref()
        .map(fs::read_to_string)
        .transpose()?;
    Ok(config)
}

fn run_call_var_bam(args: CallVarBamArgs) -> Result<()> {
    let config = call_config(&args)?;
    let reference = Reference::open(&config.reference)?;
    let output = call_var_bam(&config, args.pileup)?;
    if config.gvcf {
        fs::create_dir_all(&args.temp_file_dir)?;
        let name = format!(
            "{}.{}.{}_{}.tmp.gvcf",
            config.sample_name,
            config.contig,
            config.start.unwrap_or(1),
            config
                .end
                .unwrap_or_else(|| reference.contig_len(&config.contig) as i64)
        );
        let header = vcf_header(
            &reference,
            config.command.as_deref(),
            &config.sample_name,
            true,
            &reference.contigs()?,
        )?;
        write_gvcf(&args.temp_file_dir.join(name), &header, &output.nonvariants)?;
    }
    Ok(())
}

fn create_tensor_command(args: CallVarBamArgs, pileup: bool) -> Result<()> {
    use crate::pipeline::tensor::create_tensor_data;

    let config = call_config(&args)?;
    if args.unify_representation {
        if pileup {
            bail!("representation unification details require full-alignment input");
        }
        let reference = Reference::open(&config.reference)?;
        let details =
            crate::pipeline::unification::candidate_details_from_bam(&config, &reference)?;
        match args.unify_output {
            Some(path) if path != Path::new("PIPE") && path != Path::new("-") => {
                if let Some(parent) = path.parent() {
                    fs::create_dir_all(parent)?;
                }
                fs::write(path, details)?;
            }
            _ => print!("{details}"),
        }
        return Ok(());
    }
    let data = create_tensor_data(&config, pileup)?;
    if data.positions.is_empty() {
        return Ok(());
    }
    let reference = Reference::open(&config.reference)?;
    match args.tensor_output {
        Some(path) if path != Path::new("PIPE") && path != Path::new("-") => {
            data.write_npy(&path, &reference)?;
        }
        _ => print!("{}", data.to_text(&reference)?),
    }
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct CallVariantsArgs {
    #[arg(long = "chkpnt_fn")]
    checkpoint: PathBuf,
    #[arg(long = "ref_fn")]
    reference: PathBuf,
    #[arg(long = "sampleName", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long = "call_fn")]
    output: Option<PathBuf>,
    #[arg(long = "tensor_fn")]
    tensor: Option<PathBuf>,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long, default_value_t = false)]
    pileup: bool,
    #[arg(long, default_value_t = false)]
    add_indel_length: bool,
    #[arg(long = "showRef", action = ArgAction::SetFalse, default_value_t = true, hide = true)]
    show_reference: bool,
    #[arg(long)]
    qual: Option<f64>,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    gvcf: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    enable_long_indel: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    keep_iupac_bases: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    enable_dwell_time: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_precise: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_sensitive: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    use_gpu: bool,
    #[arg(long, default_value = "EMPTY")]
    device: String,
    #[arg(long = "cmd_fn")]
    command_file: Option<PathBuf>,
    #[arg(long = "torch_threads", default_value_t = 4)]
    torch_threads: i32,
}

fn call_variants_command(args: CallVariantsArgs) -> Result<()> {
    use crate::{
        call::decode::DecodeOptions,
        config::{MAX_LONG_VARIANT_LENGTH, MAX_VARIANT_LENGTH},
        pipeline::tensor::{TensorData, call_tensor_data},
    };

    if args.torch_threads <= 0 {
        bail!("torch_threads must be positive");
    }
    tch::set_num_threads(args.torch_threads);
    let mut data = if let Some(path) = args.tensor {
        TensorData::from_npy(&path)?
    } else {
        TensorData::from_text(
            std::io::stdin().lock(),
            args.pileup,
            args.platform,
            args.enable_dwell_time,
        )?
    };
    let reference = Reference::open(&args.reference)?;
    data.positions = data.positions_with_reference(&reference)?;
    let options = DecodeOptions {
        show_reference: args.show_reference,
        haploid_precise: args.haploid_precise,
        haploid_sensitive: args.haploid_sensitive,
        quality_for_pass: args.qual,
        add_indel_length: args.add_indel_length,
        gvcf: args.gvcf,
        pileup: args.pileup,
        enable_long_indel: args.enable_long_indel,
        maximum_variant_length: if args.enable_long_indel {
            MAX_LONG_VARIANT_LENGTH
        } else {
            MAX_VARIANT_LENGTH
        },
        keep_iupac: args.keep_iupac_bases,
    };
    let records = call_tensor_data(
        &resolve_checkpoint(&args.checkpoint)?,
        if args.pileup {
            ModelKind::Pileup
        } else {
            ModelKind::FullAlignment
        },
        &data,
        &options,
        device_from_options(args.use_gpu, &args.device)?,
    )?;
    let command = args.command_file.map(fs::read_to_string).transpose()?;
    if let Some(path) = args.output {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        crate::call::caller::write_vcf(
            &path,
            &reference,
            &records,
            &args.sample_name,
            command.as_deref(),
            args.gvcf,
        )?;
    } else {
        let header = vcf_header(
            &reference,
            command.as_deref(),
            &args.sample_name,
            args.gvcf,
            &reference.contigs()?,
        )?;
        print!("{header}");
        for record in records {
            println!("{}", record.to_line());
        }
    }
    let _ = args.contig;
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct TrainingTensorArgs {
    #[command(flatten)]
    call: CallVarBamArgs,
    #[arg(long = "var_fn")]
    truth: PathBuf,
    #[arg(long = "bin_fn")]
    output: PathBuf,
    #[arg(long, default_value_t = false)]
    allow_duplicate_chr_pos: bool,
    #[arg(long)]
    maximum_non_variant_ratio: Option<f64>,
    #[arg(long, default_value_t = false)]
    shuffle: bool,
    #[arg(long, default_value_t = false)]
    add_no_phasing_data_training: bool,
    #[arg(long, default_value_t = 0)]
    seed: u64,
}

fn create_training_tensor_command(args: TrainingTensorArgs) -> Result<()> {
    use crate::{
        pipeline::tensor::create_tensor_data,
        train::data::{TensorBinOptions, tensor_to_bin},
    };

    let mut config = call_config(&args.call)?;
    config.add_no_phasing_data_training = args.add_no_phasing_data_training;
    let data = create_tensor_data(&config, args.call.pileup)?;
    let reference = Reference::open(&config.reference)?;
    if let Some(parent) = args.output.parent() {
        fs::create_dir_all(parent)?;
    }
    let temporary = args
        .output
        .with_extension(format!("tensor.{}.tmp", std::process::id()));
    fs::write(&temporary, data.to_text(&reference)?)?;
    let result = tensor_to_bin(
        &temporary,
        &args.truth,
        args.call.bed.as_deref(),
        &args.output,
        &TensorBinOptions {
            platform: args.call.common.platform,
            pileup: args.call.pileup,
            enable_dwell_time: args.call.common.enable_dwell_time,
            allow_duplicate_positions: args.allow_duplicate_chr_pos,
            maximum_non_variant_ratio: args.maximum_non_variant_ratio,
            shuffle: args.shuffle,
            seed: args.seed,
        },
        args.call.chunk_id,
        args.call.chunk_num,
    );
    let cleanup = fs::remove_file(&temporary);
    let count = result?;
    cleanup?;
    eprintln!("[INFO] Compressed {count} tensors");
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct RealignArgs {
    #[arg(long = "bam_fn")]
    bam: PathBuf,
    #[arg(long = "ref_fn")]
    reference: PathBuf,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long = "ctgStart")]
    start: Option<i64>,
    #[arg(long = "ctgEnd")]
    end: Option<i64>,
    #[arg(long = "full_aln_regions")]
    candidate_bed: Option<PathBuf>,
    #[arg(long = "extend_bed")]
    extend_bed: Option<PathBuf>,
    #[arg(long = "call_fn")]
    output: Option<PathBuf>,
    #[arg(long = "minCoverage", default_value_t = 2.0)]
    min_coverage: f64,
    #[arg(long = "minMQ", default_value_t = 5)]
    min_mq: u8,
    #[arg(long)]
    chunk_id: Option<usize>,
    #[arg(long)]
    chunk_num: Option<usize>,
}

fn realign(args: RealignArgs) -> Result<()> {
    let reference = Reference::open(&args.reference)?;
    let output = args.output.unwrap_or_else(|| {
        args.bam.with_extension(format!(
            "{}_{}",
            args.start.unwrap_or(1),
            args.end
                .unwrap_or(reference.contig_len(&args.contig) as i64)
        ))
    });
    realign_bam(
        &args.bam,
        &reference,
        &args.contig,
        args.start,
        args.end,
        &output,
        &RealignmentOptions {
            candidate_bed: args.candidate_bed,
            extend_bed: args.extend_bed,
            min_coverage: args.min_coverage,
            min_mq: args.min_mq,
            chunk_id: args.chunk_id,
            chunk_num: args.chunk_num,
        },
    )
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct SelectQualArgs {
    #[arg(long = "output_fn")]
    output: PathBuf,
    #[arg(long, default_value_t = 0.3)]
    var_pct_full: f64,
    #[arg(long)]
    ref_pct_full: Option<f64>,
    #[arg(long)]
    var_pct_phasing: Option<f64>,
    #[arg(long, default_value_t = false)]
    phase: bool,
    #[arg(long = "vcf_fn")]
    known_vcf: Option<PathBuf>,
    #[arg(long = "qual_fn")]
    qual_name: Option<String>,
}

fn stdin_records() -> Result<Vec<VcfRecord>> {
    std::io::read_to_string(std::io::stdin())?
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(VcfRecord::from_line)
        .collect()
}

fn select_qual(args: SelectQualArgs) -> Result<()> {
    let records = stdin_records()?;
    fs::create_dir_all(&args.output)?;
    if args.phase {
        let cutoff = phase_quality(&records, args.var_pct_phasing, args.var_pct_full);
        fs::write(
            args.output
                .join(args.qual_name.unwrap_or_else(|| "phase_qual".to_string())),
            cutoff.to_string(),
        )?;
    } else {
        let (variant, reference) = full_quality(
            &records,
            args.var_pct_full,
            args.ref_pct_full,
            args.known_vcf.is_some(),
        );
        fs::write(
            args.output
                .join(args.qual_name.unwrap_or_else(|| "qual".to_string())),
            format!("{variant} {reference}"),
        )?;
    }
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct SelectHetSnpArgs {
    #[arg(long = "vcf_fn")]
    vcf: PathBuf,
    #[arg(long = "split_folder")]
    output: PathBuf,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long, default_value_t = 0.3)]
    var_pct_full: f64,
    #[arg(long = "qual_fn")]
    qual_name: Option<String>,
}

fn select_het(args: SelectHetSnpArgs) -> Result<()> {
    let file = read_vcf(&args.vcf)?;
    let qual_path = args
        .output
        .join(args.qual_name.unwrap_or_else(|| "phase_qual".to_string()));
    let cutoff = if qual_path.exists() {
        Some(fs::read_to_string(&qual_path)?.trim().parse::<f64>()?)
    } else {
        None
    };
    let selected = select_het_snp(&file.records, &args.contig, args.var_pct_full, cutoff);
    fs::create_dir_all(&args.output)?;
    let mut writer = std::io::BufWriter::new(fs::File::create(
        args.output.join(format!("{}.vcf", args.contig)),
    )?);
    writer.write_all(file.header.join("\n").as_bytes())?;
    writer.write_all(b"\n")?;
    for record in selected {
        writeln!(writer, "{}", record.to_line())?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct SelectCandidatesArgs {
    #[arg(long = "pileup_vcf_fn")]
    pileup: PathBuf,
    #[arg(long = "phased_vcf_fn")]
    phased: Option<PathBuf>,
    #[arg(long = "split_folder")]
    output: PathBuf,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long = "ref_fn")]
    reference: Option<PathBuf>,
    #[arg(long, default_value_t = 0.3)]
    var_pct_full: f64,
    #[arg(long)]
    ref_pct_full: Option<f64>,
    #[arg(long, default_value_t = 10000)]
    split_bed_size: usize,
    #[arg(long = "phasing_info_in_bam", action = ArgAction::SetFalse, default_value_t = true, hide = true)]
    phasing_info_in_bam: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    call_low_seq_entropy: bool,
    #[arg(long, default_value_t = 0.05)]
    seq_entropy_pro: f64,
    #[arg(long = "qual_fn")]
    qual_name: Option<String>,
    #[arg(long = "all_alt_fn")]
    all_alt: Option<PathBuf>,
}

fn select_full_candidates(args: SelectCandidatesArgs) -> Result<()> {
    let reference = args
        .reference
        .as_ref()
        .map(|path| Reference::open(path))
        .transpose()?;
    let qual_path = args
        .output
        .join(args.qual_name.unwrap_or_else(|| "qual".to_string()));
    let options = SelectOptions {
        var_pct_full: args.var_pct_full,
        ref_pct_full: args.ref_pct_full,
        split_bed_size: args.split_bed_size,
        phasing_info_in_bam: args.phasing_info_in_bam,
        call_low_seq_entropy: args.call_low_seq_entropy,
        seq_entropy_pro: args.seq_entropy_pro,
        qual_path: qual_path.exists().then_some(qual_path),
        all_alt_path: args.all_alt,
        ..Default::default()
    };
    select_candidates(
        &args.pileup,
        args.phased.as_deref(),
        &args.output,
        &args.contig,
        &args.platform.to_string(),
        &options,
        reference.as_ref(),
    )?;
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct MergeVcfArgs {
    #[arg(long = "pileup_vcf_fn")]
    pileup: PathBuf,
    #[arg(long = "full_alignment_vcf_fn")]
    full: PathBuf,
    #[arg(long = "output_fn")]
    output: PathBuf,
    #[arg(long = "ctgName")]
    contig: String,
    #[arg(long = "sampleName", default_value = "SAMPLE")]
    sample_name: String,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long = "bed_fn_prefix")]
    regions: Option<PathBuf>,
    #[arg(long, default_value_t = 2.0)]
    qual: f64,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    print_ref_calls: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_precise: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    haploid_sensitive: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    gvcf: bool,
    #[arg(long = "non_var_gvcf_fn")]
    nonvariants: Option<PathBuf>,
    #[arg(long = "gvcf_fn")]
    gvcf_output: Option<PathBuf>,
    #[arg(long = "ref_fn")]
    reference: Option<PathBuf>,
}

fn merge_vcf(args: MergeVcfArgs) -> Result<()> {
    let options = MergeOptions {
        print_ref_calls: args.print_ref_calls,
        haploid_precise: args.haploid_precise,
        haploid_sensitive: args.haploid_sensitive,
        qual: Some(args.qual),
    };
    let pileup = read_vcf(&args.pileup)?;
    let full = read_vcf(&args.full)?;
    let pileup: Vec<_> = pileup
        .records
        .into_iter()
        .filter(|record| record.contig == args.contig)
        .collect();
    let full: Vec<_> = full
        .records
        .into_iter()
        .filter(|record| record.contig == args.contig)
        .collect();
    let records = if args.platform == Platform::Ilmn {
        let directory = args
            .regions
            .context("Illumina merging requires --bed_fn_prefix")?;
        let regions = crate::io::interval::IntervalMap::from_bed(
            &directory.join(format!("full_aln_regions_{}", args.contig)),
            Some(&args.contig),
        )?;
        merge_illumina_records(&pileup, &full, &regions, &args.contig, &options)
    } else {
        merge_records(&pileup, &full, &options)
    };
    if let Some(parent) = args.output.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut writer = std::io::BufWriter::new(fs::File::create(&args.output)?);
    for record in &records {
        writeln!(writer, "{}", record.to_line())?;
    }
    writer.flush()?;
    if args.gvcf {
        let reference = Reference::open(
            args.reference
                .as_ref()
                .context("GVCF merging requires --ref_fn")?,
        )?;
        let nonvariants = read_vcf(
            args.nonvariants
                .as_ref()
                .context("GVCF merging requires --non_var_gvcf_fn")?,
        )?;
        let records = merge_nonvariants(
            &records,
            &nonvariants.records,
            &reference,
            args.haploid_precise || args.haploid_sensitive,
        )?;
        let header = vcf_header(
            &reference,
            None,
            &args.sample_name,
            true,
            std::slice::from_ref(&args.contig),
        )?;
        let path = args
            .gvcf_output
            .unwrap_or_else(|| args.output.with_extension("g.vcf"));
        write_gvcf(&path, &header, &records)?;
    }
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct CheckExitCodeArgs {
    #[arg(long = "parallel_log_fn")]
    log: PathBuf,
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct ModelInfoArgs {
    #[arg(long)]
    checkpoint: PathBuf,
    #[arg(long, value_enum, default_value_t = ModelKind::Pileup)]
    model_type: ModelKind,
    #[arg(long)]
    input_channels: Option<i64>,
    #[arg(long, default_value_t = false)]
    add_indel_length: bool,
    #[arg(long, default_value = "cpu")]
    device: String,
}

fn model_info(args: ModelInfoArgs) -> Result<()> {
    let model = Model::new(
        args.model_type,
        parse_device(&args.device)?,
        args.input_channels,
        args.add_indel_length,
        Some(&args.checkpoint),
    )?;
    println!(
        "{}",
        serde_json::json!({
            "checkpoint": args.checkpoint,
            "model_type": args.model_type.to_string(),
            "device": args.device,
            "tensor_count": model.var_store.len(),
        })
    );
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct TrainArgs {
    #[arg(long = "bin_fn")]
    bin_dir: PathBuf,
    #[arg(long = "ochk_prefix")]
    output_prefix: PathBuf,
    #[arg(long = "chkpnt_fn")]
    checkpoint: Option<PathBuf>,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long, default_value_t = false)]
    pileup: bool,
    #[arg(long, num_args = 0..=1, default_missing_value = "true", default_value_t = false)]
    add_indel_length: bool,
    #[arg(long, default_value_t = false)]
    enable_dwell_time: bool,
    #[arg(long, default_value_t = false)]
    reuse_bin: bool,
    #[arg(long = "maxEpoch")]
    max_epochs: Option<usize>,
    #[arg(long, default_value_t = 1e-3)]
    learning_rate: f64,
    #[arg(long)]
    exclude_training_samples: Option<String>,
    #[arg(long, default_value_t = 1)]
    mini_epochs: usize,
    #[arg(long, default_value_t = false)]
    random_validation: bool,
    #[arg(long = "validation_fn")]
    validation_dir: Option<PathBuf>,
    #[arg(long, default_value_t = -1)]
    local_rank: i32,
    #[arg(long, default_value_t = false)]
    scale_lr: bool,
    #[arg(long)]
    device: Option<String>,
    #[arg(long, hide = true)]
    seed: Option<u64>,
    #[arg(long, default_value_t = 2000, hide = true)]
    batch_size: usize,
    #[arg(long, default_value_t = 200, hide = true)]
    chunk_size: usize,
    #[arg(long = "log_fn")]
    log_path: Option<PathBuf>,
    #[arg(long)]
    best_checkpoint: Option<PathBuf>,
}

fn train_command(args: TrainArgs) -> Result<()> {
    use crate::train::trainer::{TrainingOptions, train};

    if args.pileup && args.enable_dwell_time {
        bail!("dwell time is not supported for pileup training");
    }
    if args.random_validation && args.validation_dir.is_some() {
        bail!("choose either random validation or an explicit validation directory");
    }
    let local_rank = args.local_rank >= 0;
    let distributed = local_rank || std::env::var_os("LOCAL_RANK").is_some();
    let local_rank = if local_rank {
        args.local_rank as usize
    } else if distributed {
        std::env::var("LOCAL_RANK")?.parse()?
    } else {
        0
    };
    let use_cuda = match &args.device {
        Some(device) => device.starts_with("cuda"),
        None => distributed,
    };
    let device = match &args.device {
        Some(device) => parse_device(device)?,
        None if distributed => tch::Device::Cuda(local_rank),
        None => tch::Device::cuda_if_available(),
    };
    let mut options = TrainingOptions::new(args.bin_dir, args.output_prefix, args.platform);
    options.checkpoint = args.checkpoint;
    options.pileup = args.pileup;
    options.add_indel_length = args.add_indel_length;
    options.enable_dwell_time = args.enable_dwell_time;
    options.reuse_bin = args.reuse_bin;
    options.max_epochs = args.max_epochs.unwrap_or(30);
    options.mini_epochs = args.mini_epochs;
    options.learning_rate = args.learning_rate;
    options.batch_size = args.batch_size;
    options.chunk_size = args.chunk_size;
    options.exclude_training_samples = args
        .exclude_training_samples
        .map(|values| {
            values
                .split(',')
                .filter(|value| !value.is_empty())
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default();
    options.random_validation = args.random_validation;
    options.validation_dir = args.validation_dir;
    options.device = device;
    options.seed = args.seed;
    options.scale_learning_rate = args.scale_lr;
    options.log_path = args
        .log_path
        .unwrap_or_else(|| PathBuf::from("training.log"));
    options.best_checkpoint = args
        .best_checkpoint
        .unwrap_or_else(|| PathBuf::from("best_val_loss.pt"));
    let report = if distributed {
        let world_size: usize = std::env::var("WORLD_SIZE")
            .unwrap_or_else(|_| "1".to_string())
            .parse()?;
        let rank: usize = std::env::var("RANK")
            .unwrap_or_else(|_| local_rank.to_string())
            .parse()?;
        #[cfg(feature = "distributed")]
        {
            let collective = crate::train::collective::C10dCollective::new(
                rank, world_size, local_rank, use_cuda,
            )?;
            crate::train::trainer::train_with_collective(&options, &collective)?
        }
        #[cfg(not(feature = "distributed"))]
        {
            let _ = (rank, world_size, use_cuda);
            bail!("this build does not include distributed training support");
        }
    } else {
        train(&options)?
    };
    println!(
        "{}",
        serde_json::json!({
            "epochs": report.epochs,
            "loss": report.final_loss,
            "best_validation_loss": report.best_validation_loss,
            "checkpoints": report.checkpoints,
        })
    );
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct TensorBinArgs {
    #[arg(long = "tensor_fn", default_value = "PIPE")]
    tensor: PathBuf,
    #[arg(long = "var_fn")]
    truth: PathBuf,
    #[arg(long = "bin_fn")]
    output: PathBuf,
    #[arg(long = "bed_fn")]
    bed: Option<PathBuf>,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long, default_value_t = false)]
    pileup: bool,
    #[arg(long)]
    chunk_id: Option<usize>,
    #[arg(long)]
    chunk_num: Option<usize>,
    #[arg(long, default_value_t = false)]
    allow_duplicate_chr_pos: bool,
    #[arg(long)]
    maximum_non_variant_ratio: Option<f64>,
    #[arg(long, default_value_t = false)]
    shuffle: bool,
    #[arg(long, default_value_t = false)]
    enable_dwell_time: bool,
    #[arg(long, default_value_t = 0)]
    seed: u64,
    #[arg(long = "candidate_details_fn_prefix", hide = true)]
    candidate_details: Option<String>,
}

fn tensor_bin_command(args: TensorBinArgs) -> Result<()> {
    let options = crate::train::data::TensorBinOptions {
        platform: args.platform,
        pileup: args.pileup,
        enable_dwell_time: args.enable_dwell_time,
        allow_duplicate_positions: args.allow_duplicate_chr_pos,
        maximum_non_variant_ratio: args.maximum_non_variant_ratio,
        shuffle: args.shuffle,
        seed: args.seed,
    };
    let count = crate::train::data::tensor_to_bin(
        &args.tensor,
        &args.truth,
        args.bed.as_deref(),
        &args.output,
        &options,
        args.chunk_id,
        args.chunk_num,
    )?;
    eprintln!("[INFO] Compressed {count} tensors");
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct MergeBinArgs {
    #[arg(required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,
    #[arg(long = "out_fn")]
    output: PathBuf,
    #[arg(long, default_value_t = false)]
    pileup: bool,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long, default_value_t = false)]
    enable_dwell_time: bool,
}

fn merge_bin_command(args: MergeBinArgs) -> Result<()> {
    let _ = (args.platform, args.enable_dwell_time);
    let count = crate::train::data::merge_bins(&args.inputs, &args.output, args.pileup)?;
    eprintln!("[INFO] Compressed {count} tensors");
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
struct ConvertCheckpointArgs {
    #[arg(long)]
    checkpoint: Option<PathBuf>,
    #[arg(long)]
    output: Option<PathBuf>,
    #[arg(long, value_enum, alias = "model-type")]
    model_type: Option<ModelKind>,
    #[arg(long, alias = "checkpoint-dir")]
    checkpoint_dir: Option<PathBuf>,
    #[arg(long, alias = "output-dir")]
    output_dir: Option<PathBuf>,
    #[arg(long, default_value = "pileup", alias = "pileup-prefix")]
    pileup_prefix: String,
    #[arg(long, default_value = "full_alignment", alias = "fa-prefix")]
    fa_prefix: String,
    #[arg(long, value_enum, default_value = "ont")]
    platform: Platform,
    #[arg(long, default_value_t = false, alias = "add-indel-length")]
    add_indel_length: bool,
    #[arg(long, default_value_t = false, alias = "pileup-add-indel-length")]
    pileup_add_indel_length: bool,
    #[arg(long, default_value_t = false, alias = "fa-add-indel-length")]
    fa_add_indel_length: bool,
    #[arg(long, default_value_t = false, alias = "enable-dwell-time")]
    enable_dwell_time: bool,
}

fn convert_checkpoint_command(args: ConvertCheckpointArgs) -> Result<()> {
    fn convert(
        checkpoint: PathBuf,
        output: PathBuf,
        kind: ModelKind,
        platform: Platform,
        add_indel_length: bool,
        dwell: bool,
    ) -> Result<()> {
        #[cfg(feature = "tf-checkpoint")]
        {
            crate::train::tensorflow::convert_checkpoint(
                &checkpoint,
                &output,
                kind,
                platform,
                add_indel_length,
                dwell,
            )
        }
        #[cfg(not(feature = "tf-checkpoint"))]
        {
            let _ = (checkpoint, output, kind, platform, add_indel_length, dwell);
            bail!("TensorFlow conversion requires a build with --features tf-checkpoint");
        }
    }

    if let Some(output_dir) = args.output_dir {
        let checkpoint_dir = args
            .checkpoint_dir
            .or_else(|| args.checkpoint.filter(|path| path.is_dir()))
            .context("--output_dir requires --checkpoint_dir or a checkpoint directory")?;
        fs::create_dir_all(&output_dir)?;
        convert(
            checkpoint_dir.join(&args.pileup_prefix),
            output_dir.join(format!("{}.pt", args.pileup_prefix)),
            ModelKind::Pileup,
            args.platform,
            args.pileup_add_indel_length,
            false,
        )?;
        convert(
            checkpoint_dir.join(&args.fa_prefix),
            output_dir.join(format!("{}.pt", args.fa_prefix)),
            ModelKind::FullAlignment,
            args.platform,
            args.fa_add_indel_length || args.add_indel_length,
            args.enable_dwell_time,
        )?;
        return Ok(());
    }
    convert(
        args.checkpoint.context("--checkpoint is required")?,
        args.output.context("--output is required")?,
        args.model_type.context("--model_type is required")?,
        args.platform,
        args.add_indel_length,
        args.enable_dwell_time,
    )
}

fn dispatch(cli: Cli) -> Result<()> {
    match cli.command {
        Some(Commands::Run(args)) => run_pipeline(&args.options()?),
        Some(Commands::CallVarBam(args)) => run_call_var_bam(args),
        Some(Commands::RealignReads(args)) => realign(args),
        Some(Commands::CallVariants(args)) => call_variants_command(args),
        Some(Commands::CreateTensorPileup(args)) => create_tensor_command(args, true),
        Some(Commands::CreateTensorFullAlignment(args)) => create_tensor_command(args, false),
        Some(Commands::CreateTrainingTensor(args)) => create_training_tensor_command(args),
        Some(Commands::GetTruth(args)) => crate::pipeline::training::get_truth(args),
        Some(Commands::SplitExtendBed(args)) => crate::pipeline::training::split_extend_bed(args),
        Some(Commands::UnifyRepresentation(args)) => {
            crate::pipeline::unification::unify_representation(args)
        }
        Some(Commands::SelectQual(args)) => select_qual(args),
        Some(Commands::SelectHetSnp(args)) => select_het(args),
        Some(Commands::SelectCandidates(args)) => select_full_candidates(args),
        Some(Commands::MergeVcf(args)) => merge_vcf(args),
        Some(Commands::CheckExitCode(args)) => check_parallel_log(&args.log),
        Some(Commands::CheckEnvs(args)) => crate::pipeline::commands::check_envs(args),
        Some(Commands::SortVcf(args)) => crate::pipeline::commands::sort_vcf(args),
        Some(Commands::AddBackMissingVariantsInGenotyping(args)) => {
            crate::pipeline::commands::add_back_missing(args)
        }
        Some(Commands::AddPairEndAlleleDepth(args)) => {
            crate::pipeline::commands::pair_end_depth(args)
        }
        Some(Commands::SwitchZygosityBasedOnSVCalls(args)) => {
            crate::pipeline::commands::switch_zygosity_command(args)
        }
        Some(Commands::GetOverallMetrics(args)) => {
            crate::pipeline::commands::overall_metrics_command(args)
        }
        Some(Commands::Train(args)) => train_command(args),
        Some(Commands::Tensor2Bin(args)) => tensor_bin_command(args),
        Some(Commands::MergeBin(args)) => merge_bin_command(args),
        Some(Commands::ConvertCheckpoint(args)) => convert_checkpoint_command(args),
        Some(Commands::ModelInfo(args)) => model_info(args),
        None => {
            Cli::command().print_help()?;
            println!();
            Ok(())
        }
    }
}

pub fn run() -> Result<()> {
    let mut arguments: Vec<_> = std::env::args_os().collect();
    if arguments.get(1).and_then(|argument| argument.to_str()) == Some("-v") {
        println!("Clair3 {VERSION}");
        return Ok(());
    }
    if arguments.get(1).is_some_and(|argument| {
        argument.to_str().is_some_and(|value| {
            value.starts_with('-')
                && value != "--help"
                && value != "--version"
                && value != "-h"
                && value != "-V"
        })
    }) {
        arguments.insert(1, "run".into());
    }
    dispatch(Cli::parse_from(arguments))
}
