import sys
import shlex
import subprocess
import signal
import random
import os
from os.path import dirname
from time import sleep
from argparse import ArgumentParser, SUPPRESS
import logging

logging.getLogger().setLevel(logging.INFO)


from shared.command_options import (
    CommandOption,
    CommandOptionWithNoValue,
    ExecuteCommand,
    command_string_from,
    command_option_from
)
from shared.utils import file_path_from, executable_command_string_from, subprocess_popen, str2bool, log_warning
import shared.param_p as param


class InstancesClass(object):
    def __init__(self):
        self.create_tensor = None
        self.compress_tensor = None

    def poll(self):
        self.create_tensor.poll()
        self.compress_tensor.poll()


c = InstancesClass()


def check_return_code(signum, frame):
    c.poll()
    if c.create_tensor.returncode != None and c.create_tensor.returncode != 0:
        c.compress_tensor.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if c.compress_tensor.returncode != None and c.compress_tensor.returncode != 0:
        c.create_tensor.kill()
        sys.exit("Tensor2Bin.py exited with exceptions. Exiting...")

    if (
            c.create_tensor.returncode == None or
            c.compress_tensor.returncode == None
    ):
        signal.alarm(5)


def Run(args):
    basedir = dirname(__file__)

    CTP_Bin = basedir + "/../clair3.py CreateTensorPileup"
    CTFA_Bin = basedir + "/../clair3.py CreateTensorFullAlignment"
    T2B_Bin = basedir + "/../clair3.py Tensor2Bin"

    if args.delay > 0:
        delay = random.randrange(0, args.delay)
        print("[INFO] Delay %d seconds before starting tensor creation ..." % (delay))
        sleep(delay)

    pypyBin = executable_command_string_from(args.pypy, exit_on_not_found=True)
    pythonBin = executable_command_string_from(args.python, exit_on_not_found=True)
    samtoolsBin = executable_command_string_from(args.samtools, exit_on_not_found=True)

    if args.pileup:
        bam_fn = file_path_from(args.bam_fn, exit_on_not_found=True)
    else:
        bam_fn = file_path_from(args.bam_fn)
        if bam_fn is None or bam_fn == "":
            print(log_warning(
                "[WARNING] Skip full-alignment variant calling for empty full-alignment regions"))
            return
    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    bed_fn = file_path_from(args.bed_fn)
    vcf_fn = file_path_from(args.vcf_fn)
    var_fn = file_path_from(args.var_fn, exit_on_not_found=True)
    bin_fn = args.bin_fn
    extend_bed = file_path_from(args.extend_bed)
    full_aln_regions = file_path_from(args.full_aln_regions)

    platform = args.platform
    if not platform or platform not in param.support_platform:
        sys.exit("[ERROR] Provided platform are not in support platform list [ont, hifi, ilmn]")

    pileup = args.pileup
    ctgName = args.ctgName
    min_af = args.min_af if args.min_af else param.min_af_dict[platform]
    snp_min_af = args.snp_min_af
    indel_min_af = args.indel_min_af

    if ctgName is None:
        sys.exit("--ctgName must be specified. You can call variants on multiple chromosomes simultaneously.")

    pileup_mode = command_option_from(args.pileup, 'pileup')
    phasing_info_mode = command_option_from(args.phasing_info_in_bam, 'phasing_info_in_bam')
    add_no_phasing_mode = command_option_from(args.add_no_phasing_data_training, 'add_no_phasing_data_training')
    allow_duplicate_mode = command_option_from(args.allow_duplicate_chr_pos, 'allow_duplicate_chr_pos')
    maximum_non_variant_ratio = CommandOption('maximum_non_variant_ratio', args.maximum_non_variant_ratio)
    shuffle_mode = command_option_from(args.shuffle, 'shuffle')

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

    CT_Bin = CTP_Bin if pileup else CTFA_Bin
    create_tensor_command_options = [
        pypyBin,
        CT_Bin,
        CommandOption('bam_fn', bam_fn),
        CommandOption('ref_fn', ref_fn),
        CommandOption('vcf_fn', vcf_fn),
        CommandOption('ctgName', ctgName),
        CommandOption('platform', platform),
        CommandOption('samtools', samtoolsBin),
        CommandOption('bed_fn', bed_fn),
        CommandOption('extend_bed', extend_bed),
        CommandOption('min_af', min_af),
        CommandOption('snp_min_af', snp_min_af),
        CommandOption('indel_min_af', indel_min_af),
        ctgStart,
        ctgEnd,
        chunk_id,
        chunk_num,
    ]

    if not pileup:
        create_tensor_command_options.append(phasing_info_mode)
        create_tensor_command_options.append(add_no_phasing_mode)
        create_tensor_command_options.append(CommandOption('full_aln_regions', full_aln_regions))

    compress_tensor_command_options = [
        pythonBin,
        T2B_Bin,
        CommandOption('platform', platform),
        CommandOption('var_fn', var_fn),
        CommandOption('bin_fn', bin_fn),
        CommandOption('bed_fn', bed_fn),
        chunk_id,
        chunk_num,
        allow_duplicate_mode,
        maximum_non_variant_ratio,
        shuffle_mode,
    ]
    if pileup:
        compress_tensor_command_options.append(pileup_mode)

    try:
        c.create_tensor = subprocess_popen(
            shlex.split(command_string_from(create_tensor_command_options)),
        )

        c.compress_tensor = subprocess_popen(
            shlex.split(command_string_from(compress_tensor_command_options)),
            stdin=c.create_tensor.stdout, stdout=sys.stderr
        )
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit("Failed to start required processes. Exiting...")

    signal.signal(signal.SIGALRM, check_return_code)
    signal.alarm(2)

    try:
        c.compress_tensor.wait()
        signal.alarm(0)
        c.create_tensor.stdout.close()
        c.create_tensor.wait()
    except KeyboardInterrupt as e:
        print("KeyboardInterrupt received when waiting at Tensor2Bin, terminating all scripts.")
        try:
            c.compress_tensor.terminate()
            c.create_tensor.terminate()
        except Exception as e:
            print(e)

        raise KeyboardInterrupt
    except Exception as e:
        print("Exception received when waiting at CreateTensor, terminating all scripts.")
        print(e)
        try:
            c.compress_tensor.terminate()
            c.create_tensor.terminate()
        except Exception as e:
            print(e)

        raise e


def main():
    parser = ArgumentParser(description="Create tensor binaries for pileup or full-alignment training")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="bam.bam", required=True,
                        help="BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--var_fn', type=str, default=None, required=True,
                        help="Unified VCF input filename, required")

    parser.add_argument('--bin_fn', type=str, default=None, required=True,
                        help="Compressed binary output filename, required")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--bed_fn', type=str, nargs='?', action="store", default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--min_af', type=float, default=None,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--snp_min_af', type=float, default=0.08,
                        help="Minimum SNP allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--indel_min_af', type=float, default=0.08,
                        help="Minimum Indel allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    parser.add_argument('--pypy', type=str, default="pypy3",
                        help="Path to the 'pypy', pypy3 version >= 3.6 is required, default: %(default)s")

    parser.add_argument('--python', type=str, default="python3",
                        help="Path to the 'python3', default: %(default)s")

    # options for advanced users
    parser.add_argument('--maximum_non_variant_ratio', default=None, type=float,
                        help='Maximum ratio of non-variants to variants, default: %(default)f')

    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    # options for internal process control, don't use any of them unless you are sure about the consequences
    ## In pileup mode or not
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Provide the regions to be included in full-alignment based calling
    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--add_no_phasing_data_training', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--allow_duplicate_chr_pos', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--shuffle', action='store_true',
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Wait a short while for no more than a few seconds to start the job. This is to avoid starting multiple jobs simultaneously
    ## that might use up the maximum number of threads allowed, because Tensorflow will create more threads than needed at the beginning of running the program
    ## Obseleted after adding --tensorflow_threads defaulted at 4
    parser.add_argument('--delay', type=int, default=5,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
