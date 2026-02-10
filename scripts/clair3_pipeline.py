#!/usr/bin/env python3
import argparse
import gzip
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shared.utils import str2bool


def q(value):
    return shlex.quote(str(value))


def run_command(cmd, log_path=None, env=None):
    run_env = os.environ.copy()
    run_env["PYTHONUNBUFFERED"] = "1"
    if env:
        run_env.update(env)
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=run_env,
        text=True,
        bufsize=1,
    )
    if log_path:
        log_file = open(log_path, "a", encoding="utf-8")
    else:
        log_file = None
    try:
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            if log_file:
                log_file.write(line)
                log_file.flush()
        process.wait()
    finally:
        if log_file:
            log_file.close()
    if process.returncode != 0:
        raise SystemExit(process.returncode)


def vcf_has_records(vcf_gz_path):
    try:
        with gzip.open(vcf_gz_path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line and not line.startswith("#"):
                    return True
    except FileNotFoundError:
        return False
    return False


def read_contigs(contigs_path):
    if not os.path.exists(contigs_path):
        return []
    with open(contigs_path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def add_bool_arg(parser, name, default=False):
    parser.add_argument(name, nargs="?", const=True, default=default, type=str2bool)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clair3 non-C pipeline (Python replacement for scripts/clair3.sh)"
    )
    parser.add_argument("--bam_fn", required=True)
    parser.add_argument("--ref_fn", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bed_fn", default="EMPTY")
    parser.add_argument("--vcf_fn", default="EMPTY")
    parser.add_argument("--ctg_name", default="EMPTY")
    parser.add_argument("--sample_name", default="SAMPLE")
    parser.add_argument("--chunk_num", type=int, default=0)
    parser.add_argument("--chunk_size", type=int, default=5000000)
    parser.add_argument("--qual", type=int, default=2)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--python", dest="python_bin", default="python3")
    parser.add_argument("--pypy", default="pypy3")
    parser.add_argument("--parallel", default="parallel")
    parser.add_argument("--whatshap", default="whatshap")
    parser.add_argument("--longphase", default="EMPTY")
    parser.add_argument("--var_pct_full", type=float, default=0)
    parser.add_argument("--ref_pct_full", type=float, default=0)
    parser.add_argument("--var_pct_phasing", type=float, default=0)
    parser.add_argument("--snp_min_af", type=float, default=0)
    parser.add_argument("--indel_min_af", type=float, default=0)
    parser.add_argument("--min_mq", type=int, default=5)
    parser.add_argument("--min_coverage", type=int, default=2)
    parser.add_argument("--min_contig_size", type=int, default=0)
    parser.add_argument("--pileup_model_prefix", default="pileup")
    parser.add_argument("--fa_model_prefix", default="full_alignment")
    parser.add_argument("--base_err", type=float, default=0.001)
    parser.add_argument("--gq_bin_size", type=int, default=5)

    add_bool_arg(parser, "--pileup_only", False)
    add_bool_arg(parser, "--fast_mode", False)
    add_bool_arg(parser, "--call_snp_only", False)
    add_bool_arg(parser, "--enable_variant_calling_at_sequence_head_and_tail", False)
    add_bool_arg(parser, "--output_all_contigs_in_gvcf_header", False)
    add_bool_arg(parser, "--print_ref_calls", False)
    add_bool_arg(parser, "--gvcf", False)
    add_bool_arg(parser, "--haploid_precise", False)
    add_bool_arg(parser, "--haploid_sensitive", False)
    add_bool_arg(parser, "--include_all_ctgs", False)
    add_bool_arg(parser, "--no_phasing_for_fa", False)
    add_bool_arg(parser, "--remove_intermediate_dir", False)
    add_bool_arg(parser, "--enable_phasing", False)
    add_bool_arg(parser, "--enable_long_indel", False)
    add_bool_arg(parser, "--keep_iupac_bases", False)
    add_bool_arg(parser, "--use_gpu", False)
    parser.add_argument("--device", default="EMPTY")
    add_bool_arg(parser, "--longphase_for_phasing", False)
    add_bool_arg(parser, "--use_whatshap_for_intermediate_phasing", True)
    add_bool_arg(parser, "--use_longphase_for_intermediate_phasing", False)
    add_bool_arg(parser, "--use_whatshap_for_final_output_phasing", False)
    add_bool_arg(parser, "--use_longphase_for_final_output_phasing", False)
    add_bool_arg(parser, "--use_whatshap_for_final_output_haplotagging", False)
    add_bool_arg(parser, "--enable_dwell_time", False)

    return parser.parse_args()


def main():
    args = parse_args()

    shell_folder = Path(__file__).resolve().parent
    clair3 = shell_folder / "../clair3.py"

    pileup_checkpoint = Path(args.model_path) / args.pileup_model_prefix
    if not str(pileup_checkpoint).endswith(".pt") and pileup_checkpoint.with_suffix(".pt").exists():
        pileup_checkpoint = pileup_checkpoint.with_suffix(".pt")
    full_alignment_checkpoint = Path(args.model_path) / args.fa_model_prefix
    if not str(full_alignment_checkpoint).endswith(".pt") and full_alignment_checkpoint.with_suffix(".pt").exists():
        full_alignment_checkpoint = full_alignment_checkpoint.with_suffix(".pt")

    output_folder = Path(args.output)
    log_path = output_folder / "log"
    tmp_path = output_folder / "tmp"
    split_bed_path = tmp_path / "split_beds"
    pileup_vcf_path = tmp_path / "pileup_output"
    gvcf_tmp_path = tmp_path / "gvcf_tmp_output"
    phase_output_path = tmp_path / "phase_output"
    full_alignment_output_path = tmp_path / "full_alignment_output"
    phase_vcf_path = phase_output_path / "phase_vcf"
    phase_bam_path = phase_output_path / "phase_bam"
    candidate_bed_path = full_alignment_output_path / "candidate_bed"

    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["GOTO_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"

    print("")
    print("[INFO] Check environment variables")
    cmd = (
        f"{q(args.python_bin)} {q(clair3)} CheckEnvs "
        f"--bam_fn {q(args.bam_fn)} "
        f"--bed_fn {q(args.bed_fn)} "
        f"--output_fn_prefix {q(output_folder)} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--vcf_fn {q(args.vcf_fn)} "
        f"--ctg_name {q(args.ctg_name)} "
        f"--chunk_num {q(args.chunk_num)} "
        f"--chunk_size {q(args.chunk_size)} "
        f"--include_all_ctgs {q(args.include_all_ctgs)} "
        f"--threads {q(args.threads)} "
        f"--python {q(args.python_bin)} "
        f"--pypy {q(args.pypy)} "
        f"--samtools {q(args.samtools)} "
        f"--whatshap {q(args.whatshap)} "
        f"--parallel {q(args.parallel)} "
        f"--qual {q(args.qual)} "
        f"--sampleName {q(args.sample_name)} "
        f"--var_pct_full {q(args.var_pct_full)} "
        f"--ref_pct_full {q(args.ref_pct_full)} "
        f"--snp_min_af {q(args.snp_min_af)} "
        f"--indel_min_af {q(args.indel_min_af)} "
        f"--cmd_fn {q(tmp_path / 'CMD')} "
        f"--min_contig_size {q(args.min_contig_size)} "
        f"--no_phasing_for_fa {q(args.no_phasing_for_fa)}"
    )
    run_command(cmd)

    contigs = read_contigs(tmp_path / "CONTIGS")
    if not contigs:
        print("[INFO] Exit in environment checking")
        return 0

    total_threads = args.threads

    threads_low = max(1, total_threads * 3 // 4)

    longphase_parallel = max(1, min(4, total_threads // 2))
    longphase_threads_per_job = max(1, total_threads // longphase_parallel)

    samtools_threads = max(1, total_threads // 2)

    lp_platform = "ont" if args.platform == "ont" else "pb"

    os.chdir(output_folder)

    retries = 4
    chunk_list = tmp_path / "CHUNK_LIST"

    print("[INFO] 1/7 Call variants using pileup model")
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = ""
    pileup_call_cmd = (
        f"{q(args.parallel)} --retries {retries} -C ' ' --joblog {q(log_path / 'parallel_1_call_var_bam_pileup.log')} "
        f"-j {threads_low} "
        f"\"{q(args.python_bin)} {q(clair3)} CallVarBam "
        f"--chkpnt_fn {q(pileup_checkpoint)} "
        f"--bam_fn {q(args.bam_fn)} "
        f"--call_fn {q(pileup_vcf_path / 'pileup_{1}_{2}.vcf')} "
        f"--sampleName {q(args.sample_name)} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--extend_bed {q(split_bed_path / '{1}')} "
        f"--bed_fn {q(args.bed_fn)} "
        f"--vcf_fn {q(args.vcf_fn)} "
        f"--ctgName {{1}} "
        f"--chunk_id {{2}} "
        f"--chunk_num {{3}} "
        f"--platform {q(args.platform)} "
        f"--fast_mode {q(args.fast_mode)} "
        f"--snp_min_af {q(args.snp_min_af)} "
        f"--indel_min_af {q(args.indel_min_af)} "
        f"--minMQ {q(args.min_mq)} "
        f"--minCoverage {q(args.min_coverage)} "
        f"--call_snp_only {q(args.call_snp_only)} "
        f"--enable_variant_calling_at_sequence_head_and_tail {q(args.enable_variant_calling_at_sequence_head_and_tail)} "
        f"--gvcf {q(args.gvcf)} "
        f"--base_err {q(args.base_err)} "
        f"--gq_bin_size {q(args.gq_bin_size)} "
        f"--enable_long_indel {q(args.enable_long_indel)} "
        f"--python {q(args.python_bin)} "
        f"--pypy {q(args.pypy)} "
        f"--samtools {q(args.samtools)} "
        f"--temp_file_dir {q(gvcf_tmp_path)} "
        f"--keep_iupac_bases {q(args.keep_iupac_bases)} "
        f"--cmd_fn {q(tmp_path / 'CMD')} "
        f"--pileup\" :::: {q(chunk_list)}"
    )
    run_command(pileup_call_cmd, log_path=log_path / "1_call_var_bam_pileup.log", env=env)

    run_command(
        f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_1_call_var_bam_pileup.log')}"
    )

    run_command(
        f"{q(args.pypy)} {q(clair3)} SortVcf "
        f"--input_dir {q(pileup_vcf_path)} "
        f"--vcf_fn_prefix pileup "
        f"--output_fn {q(output_folder / 'pileup.vcf')} "
        f"--sampleName {q(args.sample_name)} "
        f"--cmd_fn {q(tmp_path / 'CMD')} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--pileup_only {q(args.pileup_only)} "
        f"--print_ref_calls {q(args.print_ref_calls)} "
        f"--haploid_precise {q(args.haploid_precise)} "
        f"--haploid_sensitive {q(args.haploid_sensitive)} "
        f"--qual {q(args.qual)} "
        f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
    )

    if not vcf_has_records(output_folder / "pileup.vcf.gz"):
        print("[INFO] Exit in pileup variant calling")
        return 0

    if args.pileup_only:
        if args.remove_intermediate_dir:
            print(f"[INFO] Removing intermediate files in {tmp_path}")
            shutil.rmtree(tmp_path, ignore_errors=True)
        print(f"[INFO] Only call pileup output with --pileup_only, output file: {output_folder}/pileup.vcf.gz")
        print("[INFO] Finish calling!")
        return 0

    if args.no_phasing_for_fa:
        print("[INFO] 2/7 No phasing for full alignment calling")
        phase_bam_path.mkdir(parents=True, exist_ok=True)
        bam_path = Path(args.bam_fn)
        bam_bai = Path(f"{args.bam_fn}.bai")
        bam_bai_alt = Path(f"{Path(args.bam_fn).with_suffix('')}.bai")
        for contig in contigs:
            target_bam = phase_bam_path / f"{contig}.bam"
            if target_bam.exists() or target_bam.is_symlink():
                target_bam.unlink()
            target_bam.symlink_to(bam_path)
            if bam_bai.exists():
                target_bai = phase_bam_path / f"{contig}.bam.bai"
                if target_bai.exists() or target_bai.is_symlink():
                    target_bai.unlink()
                target_bai.symlink_to(bam_bai)
            elif bam_bai_alt.exists():
                target_bai = phase_bam_path / f"{contig}.bam.bai"
                if target_bai.exists() or target_bai.is_symlink():
                    target_bai.unlink()
                target_bai.symlink_to(bam_bai_alt)
    else:
        print("")
        print("[INFO] 2/7 Select heterozygous SNP variants for Whatshap phasing and haplotagging")
        run_command(
            f"gzip -fdc {q(output_folder / 'pileup.vcf.gz')} | {q(args.pypy)} {q(clair3)} SelectQual "
            f"--phase --output_fn {q(phase_vcf_path)} --var_pct_phasing {q(args.var_pct_phasing)}",
            log_path=log_path / "2_select_hetero_snp.log",
        )
        select_hetero_cmd = (
            f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_2_select_hetero_snp.log')} "
            f"-j {q(args.threads)} "
            f"\"{q(args.pypy)} {q(clair3)} SelectHetSnp "
            f"--vcf_fn {q(output_folder / 'pileup.vcf.gz')} "
            f"--split_folder {q(phase_vcf_path)} "
            f"--ctgName {{1}}\" ::: {' '.join(map(q, contigs))}"
        )
        run_command(select_hetero_cmd, log_path=log_path / "2_select_hetero_snp.log")
        run_command(
            f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_2_select_hetero_snp.log')}"
        )

        print("")
        if args.use_longphase_for_intermediate_phasing:
            print("[INFO] 3/7 Phase VCF file using LongPhase")
            phase_cmd = (
                f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_3_phase.log')} "
                f"-j {q(longphase_parallel)} "
                f"\"{q(args.longphase)} phase "
                f"-s {q(phase_vcf_path / '{1}.vcf')} "
                f"-b {q(args.bam_fn)} "
                f"-r {q(args.ref_fn)} "
                f"-t {q(longphase_threads_per_job)} "
                f"-o {q(phase_vcf_path / 'phased_{1}')} "
                f"--{lp_platform}\" ::: {' '.join(map(q, contigs))}"
            )
            run_command(phase_cmd, log_path=log_path / "3_phase.log")
            run_command(
                f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_3_phase.log')}"
            )
            run_command(
                f"{q(args.parallel)} -j{q(longphase_parallel)} bgzip -f {q(phase_vcf_path / 'phased_{}.vcf')} ::: {' '.join(map(q, contigs))}"
            )
        else:
            print("[INFO] 3/7 Phase VCF file using Whatshap")
            phase_cmd = (
                f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_3_phase.log')} "
                f"-j {q(args.threads)} "
                f"\"{q(args.whatshap)} phase "
                f"--output {q(phase_vcf_path / 'phased_{1}.vcf.gz')} "
                f"--reference {q(args.ref_fn)} "
                f"--chromosome {{1}} "
                f"--distrust-genotypes "
                f"--ignore-read-groups "
                f"{q(phase_vcf_path / '{1}.vcf')} "
                f"{q(args.bam_fn)}\" ::: {' '.join(map(q, contigs))}"
            )
            run_command(phase_cmd, log_path=log_path / "3_phase.log")
            run_command(
                f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_3_phase.log')}"
            )

        run_command(
            f"{q(args.parallel)} -j{q(args.threads)} tabix -f -p vcf {q(phase_vcf_path / 'phased_{}.vcf.gz')} ::: {' '.join(map(q, contigs))}"
        )

        print("")
        print("[INFO] 4/7 Haplotag input BAM file using Whatshap")
        haplotag_cmd = (
            f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_4_haplotag.log')} "
            f"-j {q(args.threads)} "
            f"\"{q(args.whatshap)} haplotag "
            f"--output {q(phase_bam_path / '{1}.bam')} "
            f"--reference {q(args.ref_fn)} "
            f"--ignore-read-groups "
            f"--regions {{1}} "
            f"{q(phase_vcf_path / 'phased_{1}.vcf.gz')} "
            f"{q(args.bam_fn)}\" ::: {' '.join(map(q, contigs))}"
        )
        run_command(haplotag_cmd, log_path=log_path / "4_haplotag.log")
        run_command(
            f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_4_haplotag.log')}"
        )
        run_command(
            f"{q(args.parallel)} -j{q(args.threads)} {q(args.samtools)} index -@{samtools_threads} {q(phase_bam_path / '{1}.bam')} ::: {' '.join(map(q, contigs))}"
        )

    print("")
    print("[INFO] 5/7 Select candidates for full-alignment calling")
    run_command(
        f"gzip -fdc {q(output_folder / 'pileup.vcf.gz')} | {q(args.pypy)} {q(clair3)} SelectQual "
        f"--output_fn {q(candidate_bed_path)} --var_pct_full {q(args.var_pct_full)} "
        f"--ref_pct_full {q(args.ref_pct_full)} --platform {q(args.platform)} --vcf_fn {q(args.vcf_fn)}",
        log_path=log_path / "5_select_candidate.log",
    )
    select_candidate_cmd = (
        f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_5_select_candidate.log')} "
        f"-j {q(args.threads)} "
        f"\"{q(args.pypy)} {q(clair3)} SelectCandidates "
        f"--pileup_vcf_fn {q(output_folder / 'pileup.vcf.gz')} "
        f"--split_folder {q(candidate_bed_path)} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--var_pct_full {q(args.var_pct_full)} "
        f"--ref_pct_full {q(args.ref_pct_full)} "
        f"--platform {q(args.platform)} "
        f"--ctgName {{1}}\" ::: {' '.join(map(q, contigs))}"
    )
    run_command(select_candidate_cmd, log_path=log_path / "5_select_candidate.log")
    run_command(
        f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_5_select_candidate.log')}"
    )

    print("")
    print("[INFO] 6/7 Call low-quality variants using full-alignment model")
    candidate_files = sorted(candidate_bed_path.glob("FULL_ALN_FILE_*") )
    if len(candidate_files) == 0:
        print("[INFO] No Candidate found! Exit in selecting full-alignment candidates")
        return 0
    full_aln_files_path = candidate_bed_path / "FULL_ALN_FILES"
    with open(full_aln_files_path, "w", encoding="utf-8") as handle:
        for path in candidate_files:
            with open(path, "r", encoding="utf-8") as infile:
                for line in infile:
                    handle.write(line.rstrip("\n") + "\n")
    full_aln_call_cmd = (
        f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_6_call_var_bam_full_alignment.log')} "
        f"-j {threads_low} "
        f"\"{q(args.python_bin)} {q(clair3)} CallVarBam "
        f"--chkpnt_fn {q(full_alignment_checkpoint)} "
        f"--bam_fn {q(phase_bam_path / '{1/.}.bam')} "
        f"--call_fn {q(full_alignment_output_path / 'full_alignment_{1/}.vcf')} "
        f"--sampleName {q(args.sample_name)} "
        f"--vcf_fn {q(args.vcf_fn)} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--full_aln_regions {{1}} "
        f"--ctgName {{1/.}} "
        f"--add_indel_length "
        f"--phasing_info_in_bam "
        f"--gvcf {q(args.gvcf)} "
        f"--minMQ {q(args.min_mq)} "
        f"--minCoverage {q(args.min_coverage)} "
        f"--enable_long_indel {q(args.enable_long_indel)} "
        f"--enable_dwell_time {q(args.enable_dwell_time)} "
        f"--python {q(args.python_bin)} "
        f"--pypy {q(args.pypy)} "
        f"--samtools {q(args.samtools)} "
        f"--keep_iupac_bases {q(args.keep_iupac_bases)} "
        f"--cmd_fn {q(tmp_path / 'CMD')} "
        f"--platform {q(args.platform)}\" :::: {q(full_aln_files_path)}"
    )
    run_command(full_aln_call_cmd, log_path=log_path / "6_call_var_bam_full_alignment.log")
    run_command(
        f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_6_call_var_bam_full_alignment.log')}"
    )

    run_command(
        f"{q(args.pypy)} {q(clair3)} SortVcf "
        f"--input_dir {q(full_alignment_output_path)} "
        f"--vcf_fn_prefix full_alignment "
        f"--output_fn {q(output_folder / 'full_alignment.vcf')} "
        f"--sampleName {q(args.sample_name)} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
    )

    if not vcf_has_records(output_folder / "full_alignment.vcf.gz"):
        print("[INFO] Exit in full-alignment variant calling")
        return 0

    if args.gvcf:
        run_command(
            f"{q(args.pypy)} {q(clair3)} SortVcf "
            f"--input_dir {q(gvcf_tmp_path)} "
            f"--vcf_fn_suffix .tmp.gvcf "
            f"--output_all_contigs_in_gvcf_header {q(args.output_all_contigs_in_gvcf_header)} "
            f"--output_fn {q(gvcf_tmp_path / 'non_var.gvcf')} "
            f"--ref_fn {q(args.ref_fn)} "
            f"--cmd_fn {q(tmp_path / 'CMD')} "
            f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
        )

    print("")
    print("[INFO] 7/7 Merge pileup VCF and full-alignment VCF")
    merge_cmd = (
        f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_7_merge_vcf.log')} "
        f"-j {q(args.threads)} "
        f"\"{q(args.pypy)} {q(clair3)} MergeVcf "
        f"--pileup_vcf_fn {q(output_folder / 'pileup.vcf.gz')} "
        f"--bed_fn_prefix {q(candidate_bed_path)} "
        f"--full_alignment_vcf_fn {q(output_folder / 'full_alignment.vcf.gz')} "
        f"--output_fn {q(tmp_path / 'merge_output/merge_{1}.vcf')} "
        f"--platform {q(args.platform)} "
        f"--print_ref_calls {q(args.print_ref_calls)} "
        f"--gvcf {q(args.gvcf)} "
        f"--haploid_precise {q(args.haploid_precise)} "
        f"--haploid_sensitive {q(args.haploid_sensitive)} "
        f"--gvcf_fn {q(tmp_path / 'merge_output/merge_{1}.gvcf')} "
        f"--non_var_gvcf_fn {q(gvcf_tmp_path / 'non_var.gvcf')} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--qual {q(args.qual)} "
        f"--ctgName {{1}}\" ::: {' '.join(map(q, contigs))}"
    )
    run_command(merge_cmd, log_path=log_path / "7_merge_vcf.log")
    run_command(
        f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_7_merge_vcf.log')}"
    )

    run_command(
        f"{q(args.pypy)} {q(clair3)} SortVcf "
        f"--input_dir {q(tmp_path / 'merge_output')} "
        f"--vcf_fn_prefix merge "
        f"--output_fn {q(output_folder / 'merge_output.vcf')} "
        f"--sampleName {q(args.sample_name)} "
        f"--ref_fn {q(args.ref_fn)} "
        f"--cmd_fn {q(tmp_path / 'CMD')} "
        f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
    )

    if not vcf_has_records(output_folder / "merge_output.vcf.gz"):
        print("[INFO] Exit in variant merging")
        return 0

    if args.gvcf:
        run_command(
            f"{q(args.pypy)} {q(clair3)} SortVcf "
            f"--input_dir {q(tmp_path / 'merge_output')} "
            f"--vcf_fn_prefix merge "
            f"--vcf_fn_suffix .gvcf "
            f"--output_all_contigs_in_gvcf_header {q(args.output_all_contigs_in_gvcf_header)} "
            f"--output_fn {q(output_folder / 'merge_output.gvcf')} "
            f"--sampleName {q(args.sample_name)} "
            f"--ref_fn {q(args.ref_fn)} "
            f"--cmd_fn {q(tmp_path / 'CMD')} "
            f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
        )

    if args.use_whatshap_for_final_output_phasing:
        print("[INFO] 7/7 Phasing VCF output in parallel using WhatsHap")
        phase_cmd = (
            f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_8_phase_vcf_output.log')} "
            f"-j {q(args.threads)} "
            f"\"{q(args.whatshap)} phase "
            f"--output {q(tmp_path / 'merge_output/phased_merge_{1}.vcf')} "
            f"--reference {q(args.ref_fn)} "
            f"--ignore-read-groups "
            f"{q(tmp_path / 'merge_output/merge_{1}.vcf')} "
            f"{q(args.bam_fn)}\" ::: {' '.join(map(q, contigs))}"
        )
        run_command(phase_cmd, log_path=log_path / "8_phase_vcf_output.log")
        run_command(
            f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_8_phase_vcf_output.log')}"
        )
        run_command(
            f"{q(args.pypy)} {q(clair3)} SortVcf "
            f"--input_dir {q(tmp_path / 'merge_output')} "
            f"--vcf_fn_prefix phased_merge "
            f"--output_fn {q(output_folder / 'phased_merge_output.vcf')} "
            f"--sampleName {q(args.sample_name)} "
            f"--ref_fn {q(args.ref_fn)} "
            f"--cmd_fn {q(tmp_path / 'CMD')} "
            f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
        )
    elif args.use_longphase_for_final_output_phasing:
        print("[INFO] 7/7 Phasing VCF output in parallel using LongPhase")
        phase_cmd = (
            f"{q(args.parallel)} --retries {retries} --joblog {q(log_path / 'parallel_8_phase_vcf_output.log')} "
            f"-j {q(longphase_parallel)} "
            f"\"{q(args.longphase)} phase "
            f"-s {q(tmp_path / 'merge_output/merge_{1}.vcf')} "
            f"-b {q(args.bam_fn)} "
            f"-r {q(args.ref_fn)} "
            f"-t {q(longphase_threads_per_job)} "
            f"-o {q(tmp_path / 'merge_output/phased_merge_{1}')} "
            f"--{lp_platform}\" ::: {' '.join(map(q, contigs))}"
        )
        run_command(phase_cmd, log_path=log_path / "3_phase.log")
        run_command(
            f"{q(args.pypy)} {q(clair3)} CheckExitCode --parallel_log_fn {q(log_path / 'parallel_8_phase_vcf_output.log')}"
        )
        run_command(
            f"{q(args.pypy)} {q(clair3)} SortVcf "
            f"--input_dir {q(tmp_path / 'merge_output')} "
            f"--vcf_fn_prefix phased_merge "
            f"--output_fn {q(output_folder / 'phased_merge_output.vcf')} "
            f"--sampleName {q(args.sample_name)} "
            f"--ref_fn {q(args.ref_fn)} "
            f"--cmd_fn {q(tmp_path / 'CMD')} "
            f"--contigs_fn {q(tmp_path / 'CONTIGS')}"
        )

    if args.use_whatshap_for_final_output_haplotagging:
        print("")
        print("[INFO] 4/7 Haplotag input BAM file using Whatshap, need some time to finish!")
        haplotag_cmd = (
            f"{q(args.whatshap)} haplotag "
            f"--output {q(output_folder / 'phased_output.bam')} "
            f"--reference {q(args.ref_fn)} "
            f"--ignore-read-groups "
            f"{q(output_folder / 'phased_merge_output.vcf.gz')} "
            f"{q(args.bam_fn)}"
        )
        run_command(haplotag_cmd, log_path=log_path / "9_haplotag.log")
        run_command(
            f"{q(args.samtools)} index -@{samtools_threads} {q(output_folder / 'phased_output.bam')}"
        )

    if args.vcf_fn != "EMPTY":
        print("[INFO] Double check re-genotyping variants")
        run_command(
            f"{q(args.pypy)} {q(clair3)} AddBackMissingVariantsInGenotyping "
            f"--vcf_fn {q(args.vcf_fn)} "
            f"--clair3_input_vcf_fn {q(output_folder / 'merge_output.vcf.gz')} "
            f"--output_fn {q(output_folder / 'merge_output.vcf')}"
        )

    if args.remove_intermediate_dir:
        print(f"[INFO] Removing intermediate files in {tmp_path}")
        shutil.rmtree(tmp_path, ignore_errors=True)

    print("")
    print(f"[INFO] Finish calling, output file: {output_folder}/merge_output.vcf.gz")
    if args.use_whatshap_for_final_output_phasing:
        print(f"[INFO] Finish calling, phased output file: {output_folder}/phased_merge_output.vcf.gz")
    if args.use_whatshap_for_final_output_haplotagging:
        print(f"[INFO] Finish calling, phased output BAM file: {output_folder}/phased_output.bam")

    return 0


if __name__ == "__main__":
    sys.exit(main())
