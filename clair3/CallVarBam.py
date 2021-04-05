import sys
import shlex
import subprocess
import multiprocessing
import signal
import random
from os.path import dirname
from time import sleep
from argparse import ArgumentParser, SUPPRESS
import logging

from shared.command_options import (
    CommandOption,
    CommandOptionWithNoValue,
    ExecuteCommand,
    command_string_from,
    command_option_from
)
from shared.utils import file_path_from, executable_command_string_from, subprocess_popen, str2bool
import shared.param_p as param


class InstancesClass(object):
    def __init__(self):
        self.create_tensor = None
        self.call_variant = None

    def poll(self):
        self.create_tensor.poll()
        self.call_variant.poll()


c = InstancesClass()


def check_return_code(signum, frame):
    c.poll()
    if c.create_tensor.returncode != None and c.create_tensor.returncode != 0:
        c.call_variant.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if c.call_variant.returncode != None and c.call_variant.returncode != 0:
        c.create_tensor.kill()
        sys.exit("call_variant.py exited with exceptions. Exiting...")

    if (
            c.create_tensor.returncode == None or
            c.call_variant.returncode == None
    ):
        signal.alarm(5)


def Run(args):
    basedir = dirname(__file__)

    CTP_Bin = basedir + "/../clair3.py CreateTensorPileup"
    CTFA_Bin = basedir + "/../clair3.py CreateTensorFullAlignment"
    RR_Bin = basedir + "/../clair3.py RealignReads"
    CVBin = basedir + "/../clair3.py CallVariants"

    pypyBin = executable_command_string_from(args.pypy, exit_on_not_found=True)
    pythonBin = executable_command_string_from(args.python, exit_on_not_found=True)
    samtoolsBin = executable_command_string_from(args.samtools, exit_on_not_found=True)

    chkpnt_fn = args.chkpnt_fn
    bam_fn = file_path_from(args.bam_fn, exit_on_not_found=True)
    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    bed_fn = file_path_from(args.bed_fn)
    extend_confident_bed_fn = file_path_from(args.extend_confident_bed_fn)
    confident_bed_fn = file_path_from(args.confident_bed_fn)

    platform = args.platform
    if not platform or platform not in param.support_platform:
        sys.exit("[ERROR] Provided platform are not in support platform list [illumina, pb, ont]")

    pileup = args.pileup
    call_fn = args.call_fn
    sampleName = args.sampleName
    ctgName = args.ctgName
    need_realignment = args.need_realignment and platform == 'illumina' and not pileup
    threshold = args.threshold if args.threshold else param.threshold_dict[platform]
    snp_min_af = args.snp_min_af
    indel_min_af = args.indel_min_af

    if ctgName is None:
        sys.exit("--ctgName must be specified. You can call variants on multiple chromosomes simultaneously.")

    haploid_precise_mode = command_option_from(args.haploid_precise, 'haploid_precise')
    haploid_sensitive_mode = command_option_from(args.haploid_sensitive, 'haploid_sensitive')
    output_for_ensemble = command_option_from(args.output_for_ensemble, 'output_for_ensemble')
    showRef_mode = command_option_from(args.showRef, 'showRef')
    qual = command_option_from(args.qual, 'qual', option_value=args.qual)

    add_indel_length_mode = CommandOption('add_indel_length', args.add_indel_length)
    phasing_info_in_bam_mode = command_option_from(args.phasing_info_in_bam, 'phasing_info_in_bam')
    need_phasing_mode = command_option_from(args.need_phasing, 'need_phasing')
    is_from_tables_mode = command_option_from(args.is_from_tables, 'is_from_tables')
    pileup_mode = command_option_from(args.pileup, 'pileup')
    gvcf_mode = CommandOption('gvcf', args.gvcf)
    fast_mode = CommandOption('fast_mode', args.fast_mode)

    ctgStart = None
    ctgEnd = None
    chunk_id = None
    chunk_num = None
    if args.ctgStart is not None and args.ctgEnd is not None and int(args.ctgStart) <= int(args.ctgEnd):
        ctgStart = CommandOption('ctgStart', args.ctgStart)
        ctgEnd = CommandOption('ctgEnd', args.ctgEnd)

    if args.chunk_id is not None and args.chunk_num is not None and int(args.chunk_id) <= int(args.chunk_num):
        chunk_id = CommandOption('chunk_id', args.chunk_id)
        chunk_num = CommandOption('chunk_num', args.chunk_num)

    if args.tensorflow_threads is None:
        numCpus = multiprocessing.cpu_count()
    else:
        numCpus = args.tensorflow_threads if args.tensorflow_threads < multiprocessing.cpu_count() else multiprocessing.cpu_count()

    maxCpus = multiprocessing.cpu_count()
    _cpuSet = ",".join(str(x) for x in random.sample(range(0, maxCpus), numCpus))

    print (_cpuSet, numCpus)
    taskSet = "taskset -c %s" % (_cpuSet)
    try:
        subprocess.check_output("which %s" % ("taskset"), shell=True)
    except:
        taskSet = ""

    if args.delay > 0:
        delay = random.randrange(0, args.delay)
        print("[INFO] Delay %d seconds before starting variant calling ..." % (delay))
        sleep(delay)

    if need_realignment:
        realign_reads_command_options = [
            pypyBin,
            RR_Bin,
            CommandOption('bam_fn', bam_fn),
            CommandOption('ref_fn', ref_fn),
            CommandOption('ctgName', ctgName),
            ctgStart,
            ctgEnd,
            chunk_id,
            chunk_num,
            CommandOption('samtools', samtoolsBin),
            CommandOption('bed_fn', bed_fn),
            CommandOption('extend_confident_bed_fn', extend_confident_bed_fn),
            CommandOption('confident_bed_fn', confident_bed_fn),
        ]
        bam_fn = "PIPE"
    CT_Bin = CTP_Bin if pileup else CTFA_Bin

    create_tensor_command_options = [
        pypyBin,
        CT_Bin,
        CommandOption('bam_fn', bam_fn),
        CommandOption('ref_fn', ref_fn),
        CommandOption('ctgName', ctgName),
        CommandOption('threshold', threshold),
        CommandOption('platform', platform),
        CommandOption('samtools', samtoolsBin),
        CommandOption('bed_fn', bed_fn),
        CommandOption('extend_confident_bed_fn', extend_confident_bed_fn),
        CommandOption('confident_bed_fn', confident_bed_fn),
        CommandOption('sampleName', args.sampleName),

        ctgStart,
        ctgEnd,
        chunk_id,
        chunk_num,
        gvcf_mode,
    ]

    if (args.gvcf):
        create_tensor_command_options.append(CommandOption('base_err', args.base_err))
        create_tensor_command_options.append(CommandOption('gq_bin_size', args.gq_bin_size))
        create_tensor_command_options.append(CommandOption('temp_file_dir', args.temp_file_dir))
        if args.bp_resolution:
            create_tensor_command_options.append(CommandOptionWithNoValue('bp_resolution'))

    if not pileup:
        create_tensor_command_options.append(phasing_info_in_bam_mode)
        create_tensor_command_options.append(need_phasing_mode)
    else:
        create_tensor_command_options.append(CommandOption('snp_min_af', snp_min_af))
        create_tensor_command_options.append(CommandOption('indel_min_af', indel_min_af))
        create_tensor_command_options.append(fast_mode)

    call_variant_command_options = [
        taskSet,
        pythonBin,
        CVBin,
        CommandOption('chkpnt_fn', chkpnt_fn),
        CommandOption('call_fn', call_fn),
        CommandOption('sampleName', sampleName),
        CommandOption('ref_fn', ref_fn),
        CommandOption('platform', platform),
        CommandOption('ctgName', ctgName),
        CommandOption('temp_file_dir', args.temp_file_dir),
        haploid_precise_mode,
        haploid_sensitive_mode,
        output_for_ensemble,
        qual,
        add_indel_length_mode,
        showRef_mode,
        is_from_tables_mode,
        pileup_mode,
        chunk_id,
        chunk_num,
        gvcf_mode,
    ]

    try:
        if need_realignment:
            c.realign_reads = subprocess_popen(
                shlex.split(command_string_from(realign_reads_command_options)),
            )
            c.create_tensor = subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options)),
                stdin=c.realign_reads.stdout)
        else:
            c.create_tensor = subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options)),
            )

        c.call_variant = subprocess_popen(
            shlex.split(command_string_from(call_variant_command_options)),
            stdin=c.create_tensor.stdout, stdout=sys.stderr
        )
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit("Failed to start required processes. Exiting...")

    signal.signal(signal.SIGALRM, check_return_code)
    signal.alarm(2)

    try:
        c.call_variant.wait()
        c.create_tensor.stdout.close()
        c.create_tensor.wait()
        if need_realignment:
            c.realign_reads.stdout.close()
            c.realign_reads.wait()
    except KeyboardInterrupt as e:
        print("KeyboardInterrupt received when waiting at CallVarBam, terminating all scripts.")
        try:
            c.call_variant.terminate()
            c.create_tensor.terminate()
            if need_realignment:
                c.realign_reads.terminate()
        except Exception as e:
            print(e)

        raise KeyboardInterrupt
    except Exception as e:
        print("Exception received when waiting at CallVarBam, terminating all scripts.")
        print(e)
        try:
            c.call_variant.terminate()
            c.create_tensor.terminate()
            if need_realignment:
                c.realign_reads.terminate()
        except Exception as e:
            print(e)

        raise e


def main():
    parser = ArgumentParser(description="Call variants using a trained model and a BAM file")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a trained model for variant calling, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, required")

    parser.add_argument('--bed_fn', type=str, nargs='?', action="store", default=None,
                        help="Call variants only in the provided bed regions, default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="bam.bam",
                        help="BAM file input, required, default: %(default)s")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="VCF output filename")

    parser.add_argument('--threshold', type=float, default=None,
                        help="Minimum allele frequence of the 1st non-reference allele for a site to be considered as a condidate site, default: %(default)f")
    
    parser.add_argument('--snp_min_af', type=float, default=0.0,
                        help="Minimum snp allele frequence for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--indel_min_af', type=float, default=0.0,
                        help="Minimum indel allele frequence for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--qual', type=int, default=None,
                        help="If set, variants with ≥$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--sampleName', type=str, nargs='?', action="store", default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required,default: %(default)s")

    parser.add_argument('--platform', type=str, default=None,
                        help="Seqeuncing platform of the input. Options: 'ont,pb,illumina', default: %(default)s")

    parser.add_argument('--pileup', action='store_true',
                        help="In pileup mode or not (full alignment mode), default: False")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--temp_file_dir', type=str, default='./',
                        help="The cache directory for storing temporary non-variant information if --gvcf is enabled, default: %(default)s")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools verision >= 1.10 is required, default: %(default)s")

    parser.add_argument('--pypy', type=str, default="pypy3",
                        help="Path to the 'pypy', pypy3 verision >= 3.6 is required, default: %(default)s")

    parser.add_argument('--python', type=str, default="python3",
                        help="Path to the 'python3', default: %(default)s")

    parser.add_argument('--delay', type=int, default=10,
                        help="Wait a short while for no more than %(default)s to start the job. This is to avoid starting multiple jobs simultaneously that might use up the maximum number of threads allowed, because Tensorflow will create more threads than needed at the beginning of running the program.")

    parser.add_argument('--showRef', action='store_false',
                        help="Show reference calls, default true for pileup model. Optional")

    parser.add_argument('--confident_bed_fn', type=str, nargs='?', action="store", default=None,
                        help="Call variant only in provied regions, works in intersection with ctgName, ctgStart and ctgEnd, optional, default: as defined by ctgName, ctgStart and ctgEnd")

    parser.add_argument('--extend_confident_bed_fn', nargs='?', action="store", type=str, default=None,
                        help="Extended regions by confident bed regions to handle mpileup with candidates near provide bed regions, default extend 16 bp distance")

    parser.add_argument('--fast_mode', type=str2bool, default=False,
                        help="Ignore low allelic frequency <= 0.15 for ont platform, default: %(default)s")

    parser.add_argument('--add_indel_length', action='store_true',
                        help="Include indel length in training and calling, false for pileup and true for raw alignment")

    # options for advanced users
    parser.add_argument('--minCoverage', type=float, default=2,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="EXPERIMENTAL: If set, reads with mapping Quality with <$minMQ will be filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="EXPERIMENTAL: If set, bases with base Quality with <$minBQ will be filtered, default: %(default)d")

    parser.add_argument('--haploid_precise', action='store_true',
                        help="EXPERIMENTAL: call haploid instead of diploid (output homo-variant only), deprecated")

    parser.add_argument('--haploid_sensitive', action='store_true',
                        help="EXPERIMENTAL: call haploid instead of diploid (output non-multi-variant only), deprecated")

    # options for debug purpose
    parser.add_argument('--ctgStart', type=int, default=None,
                        help="DEBUG: The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="DEBUG: The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Input BAM contains phasing info in HP tag, default: False")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=5, type=int,
                        help='DEBUG: Default gq bin size for merge non-variant block in gvcf option, default: %(default)d')

    parser.add_argument('--bp_resolution', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")

    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help="DEBUG: Use GPU for calling. Speed up is mostly insignficiant. Only use this for building a tested pipeline")

    parser.add_argument('--tensorflow_threads', type=int, default=4,
                        help="DEBUG: Number of threads per tensorflow job, default: %(default)s")

    # options for internal process control
    ## Output for ensemble, this mode is especially useful for extra high depth data
    parser.add_argument('--output_for_ensemble', action='store_true',
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Apply read level phasing when create tensor, will decreases calling time while costs more memoryd, default: False
    parser.add_argument('--need_phasing', action='store_true',
                        help=SUPPRESS)

    ## Apply read realignment for illumina platform, which greatly boost indel performance while increases calling time
    parser.add_argument('--need_realignment', action='store_false',
                        help=SUPPRESS)

    ## Use bin file from pytables to speed up calling, only use for testing.
    parser.add_argument('--is_from_tables', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
