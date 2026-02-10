#!/usr/bin/env python3
import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from shared.utils import str2bool

VERSION = "v2.0.0"
ERROR = "\033[31m[ERROR]\033[0m"
WARNING = "\033[33m[WARNING]\033[0m"


def add_bool_arg(parser, name, default=False):
    parser.add_argument(name, nargs="?", const=True, default=default, type=str2bool)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Clair3 entrypoint (Python replacement for run_clair3.sh). "
            "Use --bam_fn, --ref_fn, --model_path, --threads, --platform, --output."
        )
    )
    parser.add_argument("-b", "--bam_fn", required=True)
    parser.add_argument("-f", "--ref_fn", required=True)
    parser.add_argument("-m", "--model_path", required=True)
    parser.add_argument("-t", "--threads", type=int, required=True)
    parser.add_argument("-p", "--platform", required=True, choices=["ont", "hifi", "ilmn"])
    parser.add_argument("-o", "--output", required=True)

    parser.add_argument("--bed_fn", default="EMPTY")
    parser.add_argument("--vcf_fn", default="EMPTY")
    parser.add_argument("--ctg_name", default="EMPTY")
    parser.add_argument("--sample_name", default="SAMPLE")
    parser.add_argument("--qual", type=int, default=2)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--python", dest="python_bin", default="python3")
    parser.add_argument("--pypy", default="pypy3")
    parser.add_argument("--parallel", default="parallel")
    parser.add_argument("--whatshap", default="whatshap")
    parser.add_argument("--longphase", default="EMPTY")
    parser.add_argument("--chunk_num", type=int, default=0)
    parser.add_argument("--chunk_size", type=int, default=5000000)
    parser.add_argument("--snp_min_af", type=float, default=0)
    parser.add_argument("--indel_min_af", type=float, default=0)
    parser.add_argument("--var_pct_full", type=float, default=0)
    parser.add_argument("--ref_pct_full", type=float, default=0)
    parser.add_argument("--var_pct_phasing", type=float, default=0)
    parser.add_argument("--pileup_model_prefix", default="pileup")
    parser.add_argument("--fa_model_prefix", default="full_alignment")
    parser.add_argument("--min_mq", type=int, default=5)
    parser.add_argument("--min_coverage", type=int, default=2)
    parser.add_argument("--min_contig_size", type=int, default=0)
    add_bool_arg(parser, "--fast_mode", False)
    parser.add_argument("--base_err", type=float, default=0.001)
    parser.add_argument("--gq_bin_size", type=int, default=5)

    add_bool_arg(parser, "--pileup_only", False)
    add_bool_arg(parser, "--print_ref_calls", False)
    add_bool_arg(parser, "--include_all_ctgs", False)
    add_bool_arg(parser, "--gvcf", False)
    add_bool_arg(parser, "--use_whatshap_for_intermediate_phasing", True)
    add_bool_arg(parser, "--use_longphase_for_intermediate_phasing", False)
    add_bool_arg(parser, "--use_whatshap_for_final_output_phasing", False)
    add_bool_arg(parser, "--use_longphase_for_final_output_phasing", False)
    add_bool_arg(parser, "--use_whatshap_for_final_output_haplotagging", False)
    add_bool_arg(parser, "--enable_phasing", False)
    add_bool_arg(parser, "--longphase_for_phasing", False)
    add_bool_arg(parser, "--disable_c_impl", False)
    add_bool_arg(parser, "--remove_intermediate_dir", False)
    add_bool_arg(parser, "--call_snp_only", False)
    add_bool_arg(parser, "--enable_variant_calling_at_sequence_head_and_tail", False)
    add_bool_arg(parser, "--output_all_contigs_in_gvcf_header", False)
    add_bool_arg(parser, "--enable_long_indel", False)
    add_bool_arg(parser, "--keep_iupac_bases", False)
    add_bool_arg(parser, "--use_gpu", False)
    parser.add_argument("--device", default="EMPTY")
    add_bool_arg(parser, "--no_phasing_for_fa", False)
    add_bool_arg(parser, "--haploid_precise", False)
    add_bool_arg(parser, "--haploid_sensitive", False)
    add_bool_arg(parser, "--enable_dwell_time", False)

    parser.add_argument("-v", "--version", action="version", version=f"Clair3 {VERSION}")

    return parser.parse_args()


def error_exit(msg):
    print(f"{ERROR} {msg}")
    raise SystemExit(1)


def warn(msg):
    print(f"{WARNING} {msg}")


def file_exists(path):
    return Path(path).is_file()


def folder_exists(path):
    return Path(path).is_dir()


def check_bam_index(bam_path):
    return any(
        Path(p).is_file()
        for p in (
            f"{bam_path}.bai",
            f"{Path(bam_path).with_suffix('')}.bai",
            f"{bam_path}.csi",
            f"{Path(bam_path).with_suffix('')}.csi",
            f"{bam_path}.crai",
            f"{Path(bam_path).with_suffix('')}.crai",
        )
    )


def check_ref_index(ref_path):
    return any(Path(p).is_file() for p in (f"{ref_path}.fai", f"{Path(ref_path).with_suffix('')}.fai"))


def max_threads():
    detected = os.cpu_count()
    return detected if detected else 0


def ulimit_threads_limit():
    try:
        import resource

        limit = resource.getrlimit(resource.RLIMIT_NPROC)[0]
        if limit == resource.RLIM_INFINITY or limit < 0:
            return None
        return limit
    except Exception:
        return None


def model_file_exists(model_path, prefix):
    base = Path(model_path) / prefix
    return base.is_file() or base.with_suffix(".pt").is_file() or base.with_suffix(".index").is_file()


def _mv_field_has_value(field):
    parts = field.split(":", 2)
    if len(parts) != 3:
        return False

    value_type = parts[1]
    value = parts[2].strip()
    if not value:
        return False

    if value_type == "B":
        # SAM B array: TAG:B:type,val1,val2,...
        items = value.split(",", 1)
        return len(items) == 2 and bool(items[1].strip())

    return True


def check_bam_for_valid_mv_tag(samtools_path, bam_path, max_records=50):
    try:
        proc = subprocess.Popen(
            [samtools_path, "view", str(bam_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"samtools not found at '{samtools_path}'") from exc

    stderr_text = ""
    checked_records = 0
    found_mv_without_value = False

    try:
        for line in proc.stdout:
            checked_records += 1
            fields = line.rstrip("\n").split("\t")
            for field in fields[11:]:
                if not field.startswith("mv:"):
                    continue
                if _mv_field_has_value(field):
                    # Fast path: once one valid mv-tagged read is seen, stop immediately.
                    try:
                        proc.terminate()
                    except ProcessLookupError:
                        pass
                    try:
                        proc.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait()
                    return True, False, checked_records
                found_mv_without_value = True

            if checked_records >= max_records:
                # Only sample the first max_records alignments by design.
                try:
                    proc.terminate()
                except ProcessLookupError:
                    pass
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                return False, found_mv_without_value, checked_records

        # Stream ended before reaching max_records.
        stderr_text = proc.stderr.read() if proc.stderr else ""
        return_code = proc.wait()
    finally:
        if proc.stdout:
            proc.stdout.close()
        if proc.stderr:
            proc.stderr.close()

    if return_code != 0:
        detail = stderr_text.strip() or "unknown error"
        raise RuntimeError(f"samtools view failed while checking 'mv' tag: {detail}")

    return False, found_mv_without_value, checked_records


def tee_command(cmd, log_path):
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as log_file:
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_file.write(line)
            log_file.flush()
    process.wait()
    if process.returncode != 0:
        raise SystemExit(process.returncode)


def main():
    args = parse_args()

    cwd = Path.cwd()
    if str(cwd) == "/opt/bin":
        for label, value in (
            ("--bam_fn", args.bam_fn),
            ("--ref_fn", args.ref_fn),
            ("--model_path", args.model_path),
            ("--output", args.output),
        ):
            if not Path(value).is_absolute():
                error_exit(f"Require to use absolute file path {label}=FILE")
        if args.bed_fn != "EMPTY" and not Path(args.bed_fn).is_absolute():
            error_exit("Require to use absolute file path --bed_fn=FILE")
        if args.vcf_fn != "EMPTY" and not Path(args.vcf_fn).is_absolute():
            error_exit("Require to use absolute file path --vcf_fn=FILE")

    bam_path = Path(args.bam_fn)
    ref_path = Path(args.ref_fn)
    model_path = Path(args.model_path)
    output_folder = Path(args.output)

    if not bam_path.is_absolute() and bam_path.is_file():
        bam_path = bam_path.resolve()
    if not ref_path.is_absolute() and ref_path.is_file():
        ref_path = ref_path.resolve()
    if not model_path.is_absolute() and model_path.is_dir():
        model_path = model_path.resolve()
    bed_path = Path(args.bed_fn) if args.bed_fn != "EMPTY" else None
    vcf_path = Path(args.vcf_fn) if args.vcf_fn != "EMPTY" else None
    if bed_path and not bed_path.is_absolute() and bed_path.is_file():
        bed_path = bed_path.resolve()
    if vcf_path and not vcf_path.is_absolute() and vcf_path.is_file():
        vcf_path = vcf_path.resolve()
    if not output_folder.is_absolute():
        warn("No absolute output path provided, using current directory as prefix")
        output_folder = (cwd / output_folder).resolve()

    output_folder.mkdir(parents=True, exist_ok=True)
    if not output_folder.is_dir():
        error_exit(f"Cannot create output folder {output_folder}")

    ref_pro = args.ref_pct_full
    if args.platform == "ont" and ref_pro == 0:
        ref_pro = 0.1
    if args.platform != "ont" and ref_pro == 0:
        ref_pro = 0.3

    pro = args.var_pct_full
    if args.platform == "ont" and pro == 0:
        pro = 0.7
    if args.platform != "ont" and pro == 0:
        pro = 0.3

    snp_af = args.snp_min_af if args.snp_min_af != 0 else 0.08
    if args.platform == "ont" and args.indel_min_af == 0:
        indel_af = 0.15
    elif args.platform != "ont" and args.indel_min_af == 0:
        indel_af = 0.08
    else:
        indel_af = args.indel_min_af

    phasing_pct = args.var_pct_phasing if args.var_pct_phasing != 0 else 0.7
    base_model = model_path.name
    if base_model in {"r941_prom_sup_g5014", "r941_prom_hac_g5014", "ont_guppy5"}:
        phasing_pct = 0.8

    use_longphase = args.use_longphase_for_intermediate_phasing or args.longphase_for_phasing
    final_wh_phasing = args.use_whatshap_for_final_output_phasing or args.enable_phasing
    final_lp_phasing = args.use_longphase_for_final_output_phasing
    final_wh_haplotag = args.use_whatshap_for_final_output_haplotagging

    longphase_bin = args.longphase
    if (use_longphase or final_lp_phasing) and longphase_bin == "EMPTY":
        longphase_bin = str(ROOT / "longphase")
    if (use_longphase or final_lp_phasing) and not Path(longphase_bin).is_file():
        error_exit(f"Cannot find LongPhase path in {longphase_bin}, exit!")
    if (use_longphase or final_lp_phasing) and args.platform == "ilmn":
        warn("Illumina platform do not support longphase phasing, will enable whatshap phasing!")
        use_longphase = False

    if final_wh_haplotag and not final_wh_phasing and not final_lp_phasing:
        final_wh_phasing = True

    if not file_exists(bam_path):
        error_exit(f"BAM file {bam_path} not found")
    if not check_bam_index(bam_path):
        error_exit("BAM index file not found, please use 'samtools index $BAM' first")
    if not file_exists(ref_path):
        error_exit(f"Reference file {ref_path} not found")
    if not check_ref_index(ref_path):
        error_exit("Reference index fai file not found, please use 'samtools faidx $REF' first")

    if bed_path and not file_exists(bed_path):
        error_exit(f"BED file {bed_path} provides but not found")
    if vcf_path and not file_exists(vcf_path):
        error_exit(f"VCF file {vcf_path} provides but not found")
    if not folder_exists(model_path) and not os.environ.get("CONDA_PREFIX"):
        error_exit(f"Conda prefix not found, please activate clair3 conda environment first, model path: {model_path}")
    if not folder_exists(model_path):
        error_exit("Model path not found")

    max_detected = max_threads()
    if args.threads <= 0:
        error_exit("Invalid threads input --threads=INT")
    if max_detected and args.threads > max_detected:
        warn(f"Threads setting exceeds maximum available threads {max_detected}, set threads={max_detected}")
        args.threads = max_detected

    max_ulimit = ulimit_threads_limit()
    if max_ulimit:
        per_ulimit_threads = max(int(max_ulimit / 30), 1)
        if args.threads > per_ulimit_threads:
            warn(
                f"Threads setting exceeds maximum ulimit threads {args.threads} * 30 > {max_ulimit} (ulimit -u), "
                f"set threads={per_ulimit_threads}"
            )
            args.threads = per_ulimit_threads

    if args.min_mq < 5:
        warn("Invalid minimum mapping quality input --min_mq>=5")
        args.min_mq = 5
    if args.min_coverage < 2:
        warn("Invalid minimum coverage input --min_coverage>=2")
        args.min_coverage = 2
    if args.min_contig_size < 0:
        warn("Invalid minimum contig size --min_contig_size>=0")
        args.min_contig_size = 0
    if args.gq_bin_size < 0:
        warn("Invalid gq bin size --gq_bin_size>=0")
        args.gq_bin_size = 0

    if vcf_path and file_exists(vcf_path):
        snp_af = 0.0
        indel_af = 0.0

    if not model_file_exists(model_path, args.pileup_model_prefix):
        error_exit(f"No pileup model found in provided model path and model prefix {model_path}/{args.pileup_model_prefix}")
    if not model_file_exists(model_path, args.fa_model_prefix):
        error_exit(
            f"No full-alignment model found in provided model path and model prefix {model_path}/{args.fa_model_prefix}"
        )

    enable_c_impl = not args.disable_c_impl
    if enable_c_impl and args.platform == "ilmn":
        warn("Illumina platform will disable C implement to support short read realignment process!")
        enable_c_impl = False

    enable_dwell_time = args.enable_dwell_time
    if enable_dwell_time:
        if not args.platform == "ont":
            warn("Dwell time is not supported for non-ONT platform, disabling option! exit!")
            error_exit("Dwell time is not supported for non-ONT platform, disabling option! exit!")
            enable_dwell_time = False
        elif not enable_c_impl:
            warn("Dwell time requires the C implementation; re-enabling C implementation workflow.")
            enable_c_impl = True

    if enable_dwell_time:
        try:
            sampled_records = 50
            has_valid_mv, found_mv_without_value, checked_records = check_bam_for_valid_mv_tag(
                args.samtools, bam_path, max_records=sampled_records
            )
            if not has_valid_mv and found_mv_without_value:
                error_exit(
                    f"Dwell time is enabled but within the first {checked_records} alignments, "
                    "an 'mv' tag was found without a valid value. "
                    "Please ensure 'mv' is populated (e.g. mv:B:* with values) or disable --enable_dwell_time."
                )
            if not has_valid_mv:
                error_exit(
                    f"Dwell time is enabled but no valid 'mv' tag was found in the first {checked_records} alignments. "
                    "The 'mv' tag (move table) is required for dwell time analysis. "
                    "Please provide a BAM file containing 'mv' tags earlier in the file or disable --enable_dwell_time."
                )
        except FileNotFoundError:
            error_exit(f"samtools not found at '{args.samtools}', cannot verify 'mv' tag in BAM file.")
        except Exception as e:
            error_exit(f"Failed to check BAM file for 'mv' tag: {e}")

    use_gpu = args.use_gpu
    if not enable_c_impl and use_gpu:
        warn("GPU calling only support C implement for speedup!")
        use_gpu = False

    print("")
    print(f"[INFO] CLAIR3 VERSION: {VERSION}")
    print(f"[INFO] BAM FILE PATH: {bam_path}")
    print(f"[INFO] REFERENCE FILE PATH: {ref_path}")
    print(f"[INFO] MODEL PATH: {model_path}")
    print(f"[INFO] OUTPUT FOLDER: {output_folder}")
    print(f"[INFO] PLATFORM: {args.platform}")
    print(f"[INFO] THREADS: {args.threads}")
    print(f"[INFO] BED FILE PATH: {bed_path if bed_path else 'EMPTY'}")
    print(f"[INFO] VCF FILE PATH: {vcf_path if vcf_path else 'EMPTY'}")
    print(f"[INFO] CONTIGS: {args.ctg_name}")
    print(f"[INFO] CONDA PREFIX: {os.environ.get('CONDA_PREFIX', '')}")
    print(f"[INFO] SAMTOOLS PATH: {args.samtools}")
    print(f"[INFO] PYTHON PATH: {args.python_bin}")
    print(f"[INFO] PYPY PATH: {args.pypy}")
    print(f"[INFO] PARALLEL PATH: {args.parallel}")
    print(f"[INFO] WHATSHAP PATH: {args.whatshap}")
    print(f"[INFO] LONGPHASE PATH: {longphase_bin}")
    print(f"[INFO] CHUNK SIZE: {args.chunk_size}")
    if args.chunk_num > 0:
        print(f"[INFO] CHUNK NUM: {args.chunk_num}")
    if args.min_contig_size > 0:
        print(f"[INFO] MIN CONTIG SIZE: {args.min_contig_size}")
    print(f"[INFO] FULL ALIGN PROPORTION: {pro}")
    print(f"[INFO] FULL ALIGN REFERENCE PROPORTION: {ref_pro}")
    print(f"[INFO] PHASING PROPORTION: {phasing_pct}")
    print(f"[INFO] MINIMUM MQ: {args.min_mq}")
    print(f"[INFO] MINIMUM COVERAGE: {args.min_coverage}")
    print(f"[INFO] SNP AF THRESHOLD: {snp_af}")
    print(f"[INFO] INDEL AF THRESHOLD: {indel_af}")
    print(f"[INFO] BASE ERROR IN GVCF: {args.base_err}")
    print(f"[INFO] GQ BIN SIZE IN GVCF: {args.gq_bin_size}")
    print(f"[INFO] ENABLE FILEUP ONLY CALLING: {args.pileup_only}")
    print(f"[INFO] ENABLE FAST MODE CALLING: {args.fast_mode}")
    print(f"[INFO] ENABLE CALLING SNP CANDIDATES ONLY: {args.call_snp_only}")
    print(f"[INFO] ENABLE PRINTING REFERENCE CALLS: {args.print_ref_calls}")
    print(f"[INFO] ENABLE OUTPUT GVCF: {args.gvcf}")
    print(f"[INFO] ENABLE HAPLOID PRECISE MODE: {args.haploid_precise}")
    print(f"[INFO] ENABLE HAPLOID SENSITIVE MODE: {args.haploid_sensitive}")
    print(f"[INFO] ENABLE INCLUDE ALL CTGS CALLING: {args.include_all_ctgs}")
    print(f"[INFO] ENABLE NO PHASING FOR FULL ALIGNMENT: {args.no_phasing_for_fa}")
    print(f"[INFO] ENABLE REMOVING INTERMEDIATE FILES: {args.remove_intermediate_dir}")
    print(f"[INFO] ENABLE LONGPHASE FOR INTERMEDIATE VCF PHASING: {use_longphase}")
    print(f"[INFO] ENABLE PHASING FINAL VCF OUTPUT USING WHATSHAP: {final_wh_phasing}")
    print(f"[INFO] ENABLE PHASING FINAL VCF OUTPUT USING LONGPHASE: {final_lp_phasing}")
    print(f"[INFO] ENABLE HAPLOTAGGING FINAL BAM: {final_wh_haplotag}")
    print(f"[INFO] ENABLE LONG INDEL CALLING: {args.enable_long_indel}")
    print(f"[INFO] ENABLE C_IMPLEMENT: {enable_c_impl}")
    print(f"[INFO] ENABLE DWELL TIME: {enable_dwell_time}")
    if use_gpu:
        print(f"[INFO] USE GPU: {use_gpu}")
    if args.device != "EMPTY":
        print(f"[INFO] GPU DEVICE: {args.device}")

    cmd_line = " ".join(shlex.quote(item) for item in sys.argv)
    tmp_dir = output_folder / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    with open(tmp_dir / "CMD", "w", encoding="utf-8") as handle:
        handle.write(cmd_line)

    pipeline_script = "clair3_c_impl_pipeline.py" if enable_c_impl else "clair3_pipeline.py"
    pipeline_path = ROOT / "scripts" / pipeline_script

    cmd = (
        f"{shlex.quote(sys.executable)} {shlex.quote(str(pipeline_path))} "
        f"--bam_fn {shlex.quote(str(bam_path))} "
        f"--ref_fn {shlex.quote(str(ref_path))} "
        f"--threads {shlex.quote(str(args.threads))} "
        f"--model_path {shlex.quote(str(model_path))} "
        f"--platform {shlex.quote(args.platform)} "
        f"--output {shlex.quote(str(output_folder))} "
        f"--bed_fn {shlex.quote(str(bed_path) if bed_path else 'EMPTY')} "
        f"--vcf_fn {shlex.quote(str(vcf_path) if vcf_path else 'EMPTY')} "
        f"--ctg_name {shlex.quote(args.ctg_name)} "
        f"--sample_name {shlex.quote(args.sample_name)} "
        f"--chunk_num {shlex.quote(str(args.chunk_num))} "
        f"--chunk_size {shlex.quote(str(args.chunk_size))} "
        f"--samtools {shlex.quote(args.samtools)} "
        f"--python {shlex.quote(args.python_bin)} "
        f"--pypy {shlex.quote(args.pypy)} "
        f"--parallel {shlex.quote(args.parallel)} "
        f"--whatshap {shlex.quote(args.whatshap)} "
        f"--qual {shlex.quote(str(args.qual))} "
        f"--var_pct_full {shlex.quote(str(pro))} "
        f"--ref_pct_full {shlex.quote(str(ref_pro))} "
        f"--var_pct_phasing {shlex.quote(str(phasing_pct))} "
        f"--snp_min_af {shlex.quote(str(snp_af))} "
        f"--indel_min_af {shlex.quote(str(indel_af))} "
        f"--min_mq {shlex.quote(str(args.min_mq))} "
        f"--min_coverage {shlex.quote(str(args.min_coverage))} "
        f"--min_contig_size {shlex.quote(str(args.min_contig_size))} "
        f"--pileup_only {shlex.quote(str(args.pileup_only))} "
        f"--gvcf {shlex.quote(str(args.gvcf))} "
        f"--base_err {shlex.quote(str(args.base_err))} "
        f"--gq_bin_size {shlex.quote(str(args.gq_bin_size))} "
        f"--fast_mode {shlex.quote(str(args.fast_mode))} "
        f"--call_snp_only {shlex.quote(str(args.call_snp_only))} "
        f"--enable_variant_calling_at_sequence_head_and_tail {shlex.quote(str(args.enable_variant_calling_at_sequence_head_and_tail))} "
        f"--output_all_contigs_in_gvcf_header {shlex.quote(str(args.output_all_contigs_in_gvcf_header))} "
        f"--print_ref_calls {shlex.quote(str(args.print_ref_calls))} "
        f"--haploid_precise {shlex.quote(str(args.haploid_precise))} "
        f"--haploid_sensitive {shlex.quote(str(args.haploid_sensitive))} "
        f"--include_all_ctgs {shlex.quote(str(args.include_all_ctgs))} "
        f"--no_phasing_for_fa {shlex.quote(str(args.no_phasing_for_fa))} "
        f"--pileup_model_prefix {shlex.quote(args.pileup_model_prefix)} "
        f"--fa_model_prefix {shlex.quote(args.fa_model_prefix)} "
        f"--remove_intermediate_dir {shlex.quote(str(args.remove_intermediate_dir))} "
        f"--enable_phasing {shlex.quote(str(final_wh_phasing))} "
        f"--enable_long_indel {shlex.quote(str(args.enable_long_indel))} "
        f"--keep_iupac_bases {shlex.quote(str(args.keep_iupac_bases))} "
        f"--use_gpu {shlex.quote(str(use_gpu))} "
        f"--device {shlex.quote(args.device)} "
        f"--longphase_for_phasing {shlex.quote(str(use_longphase))} "
        f"--longphase {shlex.quote(longphase_bin)} "
        f"--use_whatshap_for_intermediate_phasing {shlex.quote(str(args.use_whatshap_for_intermediate_phasing))} "
        f"--use_longphase_for_intermediate_phasing {shlex.quote(str(use_longphase))} "
        f"--use_whatshap_for_final_output_phasing {shlex.quote(str(final_wh_phasing))} "
        f"--use_longphase_for_final_output_phasing {shlex.quote(str(final_lp_phasing))} "
        f"--use_whatshap_for_final_output_haplotagging {shlex.quote(str(final_wh_haplotag))} "
        f"--enable_dwell_time {shlex.quote(str(enable_dwell_time))}"
    )

    log_path = output_folder / "run_clair3.log"
    tee_command(cmd, log_path)


if __name__ == "__main__":
    main()
