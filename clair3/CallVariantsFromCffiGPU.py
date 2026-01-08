import sys
import logging
from argparse import ArgumentParser, SUPPRESS
import subprocess
import os

from shared.utils import str2bool, log_error, log_warning, str_none
logging.basicConfig(format='%(message)s', level=logging.INFO)

file_directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
main_entry = os.path.join(file_directory, "clair3.py")

def get_gpu_memory(gpu_id):
    command = "nvidia-smi --query-gpu=memory.free --format=csv "
    if gpu_id is not None:
        command += f" --id={gpu_id}"
    memory_free_info = subprocess.check_output(command.split()).decode('ascii').split('\n')[:-1][1:]
    memory_free_values = [int(x.split()[0]) for i, x in enumerate(memory_free_info)]
    return memory_free_values

def check_gpu_memory(memory, device_ids=None, print_log=True):
    import torch
    all_device_ids = list(range(torch.cuda.device_count()))
    if device_ids is None:
        device_ids = all_device_ids
    if not all_device_ids:
        return
    gpu_id_list = []
    detect_gpu = False
    for device_id in all_device_ids:
        detect_gpu = True
        free_mem = get_gpu_memory(gpu_id=device_id)[0]  # Convert to MB
        gpu_threads = free_mem // memory
        gpu_id_list += [device_id] * gpu_threads
        if print_log:
            print(f"GPU {device_id} free memory: {free_mem} MB, assigning {memory} MB per thread, {gpu_threads} threads available")
    if len(device_ids) == 0 or not detect_gpu:
        print(log_error("No GPU available, Please disabling --use_gpu for variant calling, exiting."))
        sys.exit(1)
    if len(gpu_id_list) == 0:
        print(log_error("No memory in GPU, Please assign GPU memory first, exiting."))
        sys.exit(1)
    return gpu_id_list

def Run(args):
    try:
        rc = subprocess.check_call('time', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        time = 'time '
    except subprocess.CalledProcessError as e:
        time = ''

    prefix = "full_alignment" if not args.pileup else "pileup"
    dwell_flag = ' --enable_dwell_time True' if args.enable_dwell_time else ''

    pileup_per_thread_gpu_memory = 5000  # GB
    full_alignment_per_thread_gpu_memory = 8000  # GB

    cuda_visual_devices = os.environ.get("CUDA_VISIBLE_DEVICES", "-1")

    if args.device and args.device != "EMPTY":
        if not args.device.startswith('cuda:'):
            print(log_error("Please specify the device as '--device=cuda:0' or '--device=cuda:0,1'"))
            sys.exit(1)
        device_ids = [int(x) for x in args.device.split("cuda:")[-1].split(",")]
        os.environ["CUDA_VISIBLE_DEVICES"] = ','.join([str(item) for item in device_ids])
    else:
        device_ids = None
        #check if there is any preset cuda env
        if cuda_visual_devices != "-1":
            cuda_device_ids = [int(x) for x in cuda_visual_devices.split(",") if x.isdigit()]
            print(f"[INFO] ENV 'CUDA_VISIBLE_DEVICES' is set to {cuda_device_ids}, using these GPUs only.")

    gpu_id_list = check_gpu_memory(pileup_per_thread_gpu_memory if args.pileup else full_alignment_per_thread_gpu_memory, device_ids, print_log=False)

    if args.pileup:
        cp_command = time + args.parallel + ' -C " " '
        cp_command += ' --joblog ' + args.output_dir + '/log/parallel_1_pileup_create_tensor.log'
        cp_command += ' -j ' + str(args.cpu_threads)
        cp_command += ' --retries 4'
        cp_command += ' ' + args.python + ' ' + main_entry + " CallVariantsFromCffi"
        cp_command += ' --chkpnt_fn ' + args.chkpnt_fn
        cp_command += ' --bam_fn ' + args.bam_fn
        cp_command += " --call_fn " + os.path.join(args.call_fn, "pileup_{1}_{2}_call.vcf.gz")
        cp_command += ' --sampleName ' + args.sampleName
        cp_command += ' --ref_fn ' + args.ref_fn if args.ref_fn else ''
        cp_command += ' --vcf_fn ' + args.vcf_fn if args.vcf_fn != "EMPTY" else ''
        cp_command += ' --extend_bed ' + args.output_dir + "/tmp/split_beds/{1}"
        if args.bed_fn != "" and args.bed_fn is not None:
            cp_command += ' --bed_fn ' + str(args.bed_fn)
        cp_command += ' --threads ' + str(args.internal_threads)
        cp_command += ' --ctgName {1} '
        cp_command += ' --chunk_id {2} '
        cp_command += ' --chunk_num {3} '
        cp_command += ' --platform ' + args.platform
        cp_command += ' --fast_mode ' + str(args.fast_mode)
        cp_command += ' --snp_min_af ' + str(args.snp_min_af)
        cp_command += ' --indel_min_af ' + str(args.indel_min_af)
        cp_command += ' --minMQ ' + str(args.minMQ)
        cp_command += ' --minCoverage ' + str(args.minCoverage)
        cp_command += ' --call_snp_only ' + str(args.call_snp_only)
        cp_command += ' --enable_variant_calling_at_sequence_head_and_tail ' + str(
            args.enable_variant_calling_at_sequence_head_and_tail)
        cp_command += ' --gvcf ' + str(args.gvcf)
        cp_command += ' --base_err ' + str(args.base_err)
        cp_command += ' --gq_bin_size ' + str(args.gq_bin_size)
        cp_command += ' --enable_long_indel ' + str(args.enable_long_indel)
        cp_command += ' --samtools ' + str(args.samtools)
        cp_command += ' --temp_file_dir ' + args.output_dir + "/tmp/gvcf_tmp_output"
        cp_command += ' --pileup '
        cp_command += ' --keep_iupac_bases ' + str(args.keep_iupac_bases)
        cp_command += ' --cmd_fn ' + str(args.cmd_fn)
        cp_command += ' --use_gpu ' + str(args.use_gpu)
        cp_command += " --tensor_can_fn " + args.output_dir + f"/tmp/pileup_output/{prefix}" + "_{1}_{2}"
        cp_command += ' :::: ' + args.output_dir + "/tmp/CHUNK_LIST"

        try:
            return_code = subprocess.check_call(cp_command, shell=True, stdout=sys.stdout)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                "ERROR in Pileup tensor generation, THE FOLLOWING COMMAND FAILED: {}\n".format(cp_command))
            exit(1)

        file_list = []
        ctg_set = set()
        with open(args.output_dir + "/tmp/CHUNK_LIST") as f:
            for row in f:
                ctg, chunk_id, chunk_num = row.rstrip().split(' ')
                if os.path.exists(
                        os.path.join(args.output_dir, 'tmp', 'pileup_output', f'{prefix}_{ctg}_{chunk_id}.npy')):
                    file_list.append(
                        os.path.join(args.output_dir, 'tmp', 'pileup_output', f'{prefix}_{ctg}_{chunk_id}'))
                    ctg_set.add(ctg)
        if len(file_list) == 0:
            return
        logging.info("[INFO] Found {} pileup tensor files to process".format(len(file_list)))


        # check gpu resource again, split the file based on the GPU resource
        gpu_id_list = check_gpu_memory(pileup_per_thread_gpu_memory, device_ids)

        total_gpu_treads = len(gpu_id_list)
        each_gpu_file_count = len(file_list) // total_gpu_treads if len(file_list) % total_gpu_treads == 0 else len(
            file_list) // total_gpu_treads + 1

        all_gpu_chunk_file_list = os.path.join(args.output_dir, 'tmp', f'{prefix}_all_gpu_chunk')
        with open(all_gpu_chunk_file_list, 'w') as all_f:
            for i in range(total_gpu_treads):
                gpu_file_list = file_list[i * each_gpu_file_count: (i + 1) * each_gpu_file_count]
                if len(gpu_file_list) == 0:
                    continue

                gpu_chunk_num_fn = os.path.join(args.output_dir, 'tmp', f'{prefix}_gpu_chunk_{i}')
                with open(gpu_chunk_num_fn, 'w') as f:
                    for fn in gpu_file_list:
                        f.write(fn + '\n')
                all_f.write(' '.join([str(i), str(gpu_id_list[i]), gpu_chunk_num_fn]) + '\n')


        for f in os.listdir(args.call_fn):
            if f.startswith(prefix) and f.endswith('.vcf'):
                os.remove(os.path.join(args.call_fn, f))

        cp_command = time + args.parallel + ' -C " " '
        cp_command += ' --joblog ' + args.output_dir + '/log/parallel_1_pileup_call_variant.log'
        cp_command += ' -j ' + str(total_gpu_treads)
        cp_command += ' --retries 4'
        cp_command += ' ' + args.python + ' ' + main_entry + " CallVariantsFromCffi"
        cp_command += ' --chkpnt_fn ' + args.chkpnt_fn
        cp_command += ' --bam_fn ' + args.bam_fn
        cp_command += " --call_fn " + os.path.join(args.call_fn, f"{prefix}_" + "{1}.vcf")
        cp_command += ' --sampleName ' + args.sampleName
        cp_command += ' --ref_fn ' + args.ref_fn if args.ref_fn else ''
        cp_command += ' --vcf_fn ' + args.vcf_fn if args.vcf_fn != "EMPTY" else ''
        if args.bed_fn != "" and args.bed_fn is not None:
            cp_command += ' --bed_fn ' + str(args.bed_fn)
        cp_command += ' --threads ' + str(args.internal_threads)
        cp_command += ' --platform ' + args.platform
        cp_command += ' --fast_mode ' + str(args.fast_mode)
        cp_command += ' --snp_min_af ' + str(args.snp_min_af)
        cp_command += ' --indel_min_af ' + str(args.indel_min_af)
        cp_command += ' --minMQ ' + str(args.minMQ)
        cp_command += ' --minCoverage ' + str(args.minCoverage)
        cp_command += ' --call_snp_only ' + str(args.call_snp_only)
        cp_command += ' --enable_variant_calling_at_sequence_head_and_tail ' + str(
            args.enable_variant_calling_at_sequence_head_and_tail)
        cp_command += ' --gvcf ' + str(args.gvcf)
        cp_command += ' --base_err ' + str(args.base_err)
        cp_command += ' --gq_bin_size ' + str(args.gq_bin_size)
        cp_command += ' --enable_long_indel ' + str(args.enable_long_indel)
        cp_command += ' --samtools ' + str(args.samtools)
        cp_command += ' --temp_file_dir ' + args.output_dir + "/tmp/gvcf_tmp_output"
        cp_command += ' --pileup '
        cp_command += ' --keep_iupac_bases ' + str(args.keep_iupac_bases)
        cp_command += ' --cmd_fn ' + str(args.cmd_fn)
        cp_command += ' --gpu_id {2} '
        cp_command += ' --cpu_threads ' + str(args.cpu_threads)
        cp_command += ' --use_gpu ' + str(args.use_gpu)
        cp_command += " --output_tensor_can_fn_list {3} "
        cp_command += ' :::: ' + all_gpu_chunk_file_list

        try:
            return_code = subprocess.check_call(cp_command, shell=True, stdout=sys.stdout)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                "ERROR in Pileup model calling, THE FOLLOWING COMMAND FAILED: {}\n".format(cp_command))
            exit(1)

        print("[INFO] Removing temporary tensor files...")
        for f in file_list:
            if os.path.exists(f + ".npy"):
                os.remove(f + ".npy")
            if os.path.exists(f + ".info"):
                os.remove(f + ".info")

    else:
        # create full-alignment tensor
        ct_command = '(' + time + args.parallel + ' -C " "'
        ct_command += ' --joblog ' + args.output_dir + '/log/parallel_6_full_alignment_create_tensor.log'
        ct_command += ' -j ' + str(args.cpu_threads)
        ct_command += ' --retries 4'
        ct_command += ' ' + args.python + ' ' + main_entry + " CallVariantsFromCffi"
        ct_command += ' --chkpnt_fn ' + args.chkpnt_fn
        ct_command += ' --bam_fn ' + args.bam_fn
        ct_command += ' --sampleName ' + args.sampleName
        ct_command += ' --vcf_fn ' + args.vcf_fn if args.vcf_fn != "EMPTY" else ''
        ct_command += ' --ref_fn ' + args.ref_fn if args.ref_fn else ''
        ct_command += ' --full_aln_regions {1}'
        ct_command += ' --ctgName {1/.} '
        ct_command += " --add_indel_length" if args.add_indel_length else ''
        ct_command += ' --no_phasing_for_fa ' + str(args.no_phasing_for_fa) if args.no_phasing_for_fa else ''
        ct_command += ' --minMQ ' + str(args.minMQ)
        ct_command += ' --minCoverage ' + str(args.minCoverage)
        ct_command += ' --phased_vcf_fn ' + os.path.join(args.phased_vcf_fn,
                                                         "phased_{1/.}.vcf.gz") if args.phased_vcf_fn else ''
        ct_command += ' --gvcf ' + str(args.gvcf)
        ct_command += ' --enable_long_indel ' + str(args.enable_long_indel)
        ct_command += ' --samtools ' + args.samtools
        ct_command += ' --use_gpu ' + str(args.use_gpu)
        ct_command += ' --keep_iupac_bases ' + str(args.keep_iupac_bases)
        ct_command += ' --cmd_fn ' + args.cmd_fn if args.cmd_fn else ''
        ct_command += ' --platform ' + args.platform
        ct_command += ' --use_gpu ' + str(args.use_gpu)
        ct_command += dwell_flag
        ct_command += " --tensor_can_fn " + args.output_dir + "/tmp/full_alignment_output/{1/}"
        ct_command += ' :::: ' + args.full_aln_files
        ct_command += ' ) 2>&1 | tee ' + args.output_dir + '/log/6_full_alignment_create_tensor.log'

        try:
            return_code = subprocess.check_call(ct_command, shell=True, stdout=sys.stdout)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                "ERROR in Full-alignment create tensor, THE FOLLOWING COMMAND FAILED: {}\n".format(ct_command))
            exit(1)

        # check gpu resource again, split the file based on the GPU resource
        gpu_id_list = check_gpu_memory(full_alignment_per_thread_gpu_memory, device_ids)


        file_list = []
        for f in open(args.full_aln_files, 'r').readlines():
            f = f.strip()
            if not f:
                continue
            suffix = f.split('/')[-1]
            if os.path.exists(args.output_dir + f"/tmp/full_alignment_output/{suffix}.npy"):
                file_list.append(args.output_dir + f"/tmp/full_alignment_output/{suffix}")

        total_gpu_treads = len(gpu_id_list)
        each_gpu_file_count = len(file_list) // total_gpu_treads if len(file_list) % total_gpu_treads == 0 else len(
            file_list) // total_gpu_treads + 1

        all_gpu_chunk_file_list = os.path.join(args.output_dir, 'tmp', f'{prefix}_all_gpu_chunk')
        with open(all_gpu_chunk_file_list, 'w') as all_f:
            for i in range(total_gpu_treads):
                gpu_file_list = file_list[i * each_gpu_file_count: (i + 1) * each_gpu_file_count]
                if len(gpu_file_list) == 0:
                    continue

                gpu_chunk_num_fn = os.path.join(args.output_dir, 'tmp', f'{prefix}_gpu_chunk_{i}')
                with open(gpu_chunk_num_fn, 'w') as f:
                    for fn in gpu_file_list:
                        f.write(fn + '\n')
                all_f.write(' '.join([str(i), str(gpu_id_list[i]), gpu_chunk_num_fn]) + '\n')

        for f in os.listdir(args.call_fn):
            if f.startswith(prefix) and f.endswith('.vcf'):
                os.remove(os.path.join(args.call_fn, f))

        #full-alignment variant calling
        ct_command = '(' + time + args.parallel + ' -C " "'
        ct_command += ' --joblog ' + args.output_dir + '/log/parallel_6_full_alignment_call_variant.log'
        ct_command += ' -j ' + str(total_gpu_treads)
        ct_command += ' --retries 4'
        ct_command += ' ' + args.python + ' ' + main_entry + " CallVariantsFromCffi"
        ct_command += ' --chkpnt_fn ' + args.chkpnt_fn
        ct_command += ' --bam_fn ' + args.bam_fn
        ct_command += ' --call_fn ' + os.path.join(args.call_fn, f"{prefix}_" + "{1}.vcf")
        ct_command += ' --sampleName ' + args.sampleName
        ct_command += ' --vcf_fn ' + args.vcf_fn if args.vcf_fn != "EMPTY" else ''
        ct_command += ' --ref_fn ' + args.ref_fn if args.ref_fn else ''
        ct_command += " --add_indel_length" if args.add_indel_length else ''
        ct_command += ' --no_phasing_for_fa ' + str(args.no_phasing_for_fa) if args.no_phasing_for_fa else ''
        ct_command += ' --minMQ ' + str(args.minMQ)
        ct_command += ' --minCoverage ' + str(args.minCoverage)
        ct_command += ' --phased_vcf_fn ' + os.path.join(args.phased_vcf_fn,
                                                         "phased_{2}.vcf.gz") if args.phased_vcf_fn else ''
        ct_command += ' --gvcf ' + str(args.gvcf)
        ct_command += ' --enable_long_indel ' + str(args.enable_long_indel)
        ct_command += ' --samtools ' + args.samtools
        ct_command += ' --use_gpu ' + str(args.use_gpu)
        ct_command += ' --keep_iupac_bases ' + str(args.keep_iupac_bases)
        ct_command += ' --cmd_fn ' + args.cmd_fn if args.cmd_fn else ''
        ct_command += ' --platform ' + args.platform
        ct_command += ' --gpu_id {2}'
        ct_command += ' --cpu_threads ' + str(args.cpu_threads)
        ct_command += ' --use_gpu ' + str(args.use_gpu)
        ct_command += dwell_flag
        ct_command += ' --output_tensor_can_fn_list {3} '
        ct_command += ' :::: ' + all_gpu_chunk_file_list
        ct_command += ' ) 2>&1 | tee ' + args.output_dir + '/log/6_full_alignment_call_variant.log'

        try:
            return_code = subprocess.check_call(ct_command, shell=True, stdout=sys.stdout)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(
                "ERROR in Full-alignment model calling, THE FOLLOWING COMMAND FAILED: {}\n".format(ct_command))
            exit(1)

        print("[INFO] Removing temporary tensor files...")
        for f in file_list:
            if os.path.exists(f + ".npy"):
                os.remove(f + ".npy")
            if os.path.exists(f + ".info"):
                os.remove(f + ".info")

def main():
    parser = ArgumentParser(description="The wrapper of GPU Call variants")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--full_aln_files', type=str, default=None,
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Tensor input filename, or stdin if not set")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a trained model for variant calling, required")

    parser.add_argument('--call_fn', type=str, default="clair3",
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required if --gvcf is enabled")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--qual', type=int, default=2,
                        help="If set, variants with >=$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    # options for advanced users
    parser.add_argument('--temp_file_dir', type=str, default='./',
                        help="EXPERIMENTAL: The cache directory for storing temporary non-variant information if --gvcf is enabled, default: %(default)s")

    parser.add_argument('--haploid_precise', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant")

    parser.add_argument('--haploid_sensitive', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant")

    parser.add_argument('--enable_long_indel', type=str2bool, default=False,
                        help="EXPERIMENTAL: Enable long Indel variants(>50 bp) calling")

    parser.add_argument('--keep_iupac_bases', type=str2bool, default=False,
                        help="EXPERIMENTAL: Keep IUPAC (non ACGTN) reference and alternate bases, default: convert all IUPAC bases to N")

    parser.add_argument('--enable_dwell_time', type=str2bool, default=False,
                        help="EXPERIMENTAL: Enable dwell time channel for full-alignment tensors")

    # options for debug purpose
    parser.add_argument('--use_gpu', type=str2bool, default=True,
                        help="Use GPU for calling")

    parser.add_argument('--threads', type=int, default=1,
                        help="How many threads for python to use for parallelization default: %(default)f")

    parser.add_argument('--predict_fn', type=str, default=None,
                        help="DEBUG: Output network output probabilities for further analysis")

    parser.add_argument('--input_probabilities', action='store_true',
                        help="DEBUG: Use network probability outputs as input and generate variants from them")

    parser.add_argument('--output_probabilities', action='store_true',
                        help="DEBUG: Output the network probabilities of gt21, genotype, indel_length_1 and indel_length_2")

    parser.add_argument('--enable_variant_calling_at_sequence_head_and_tail', type=str2bool, default=False,
                        help="EXPERIMENTAL: Enable variant calling in sequence head and tail start or end regions that flanking 16bp windows having no read support. Default: disable.")

    # options for internal process control
    ## In pileup mode or not (full alignment mode), default: False
    parser.add_argument('--pileup', action='store_true',
                        help=SUPPRESS)

    ## Include indel length in training and calling, false for pileup and true for raw alignment
    parser.add_argument('--add_indel_length', action='store_true',
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Enable debug mode, default: False, optional
    parser.add_argument('--debug', action='store_true',
                        help=SUPPRESS)

    ## Generating outputs for ensemble model calling
    parser.add_argument('--output_for_ensemble', action='store_true',
                        help=SUPPRESS)

    ## Use bin file from HDF5 to speed up calling.
    parser.add_argument('--is_from_tables', action='store_true',
                        help=SUPPRESS)

    ## Output reference calls
    parser.add_argument('--showRef', type=str2bool, default=True,
                        help=SUPPRESS)

    # Pileup create tensor options for pileup calling
    parser.add_argument('--bam_fn', type=str, default="input.bam", required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--bed_fn', type=str, nargs='?', action="store", default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--snp_min_af', type=float, default=0.08,
                        help="Minimum snp allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--indel_min_af', type=float, default=0.15,
                        help="Minimum indel allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--minCoverage', type=int, default=2,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$minMQ are filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="EXPERIMENTAL: If set, bases with base quality with <$minBQ are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=144,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    parser.add_argument('--fast_mode', type=str2bool, default=False,
                        help="EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: %(default)s")

    parser.add_argument('--call_snp_only', type=str2bool, default=False,
                        help="EXPERIMENTAL: Call candidates pass snp minimum AF only, ignore Indel candidates")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=5, type=int,
                        help='DEBUG: Default gq bin size for merge non-variant block in gvcf option, default: %(default)d')

    parser.add_argument('--bp_resolution', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")

    # Full-alignment create tensor options for full-alignment calling
    parser.add_argument('--phased_vcf_fn', type=str, default=None,
                        help="Use heterozygous SNP variants in phased vcf file for haplotaging")

    parser.add_argument('--no_phasing_for_fa', type=str2bool, default=False,
                        help="EXPERIMENTAL: Call variants without whatshap or longphase phasing in full alignment calling")

    ## Provide the regions to be included in full-alignment based calling
    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help=SUPPRESS)

    ## If defined, added command line into VCF header
    parser.add_argument('--cmd_fn', type=str_none, default=None,
                        help=SUPPRESS)

    parser.add_argument(
        "--parallel", type=str, default="parallel",
                        help="Absolute path of parallel, parallel >= 20191122 is required."
    )

    parser.add_argument(
        "--pypy", type=str, default="pypy3",
                        help="Absolute path of pypy3, pypy3 >= 3.6 is required."
    )

    parser.add_argument(
        "--python", type=str, default="python3",
                        help="Absolute path of pypy3, pypy3 >= 3.6 is required."
    )

    parser.add_argument(
        '--cpu_threads', type=int, default=4,
                        help="How many threads for python to use for parallelization default: %(default)f")

    parser.add_argument(
        '--gpu_threads', type=int, default=4,
                        help="How many threads for python to use for parallelization default: %(default)f")

    parser.add_argument(
        '--internal_threads', type=int, default=1,
                        help="How many threads for python to use for parallelization default: %(default)f")

    parser.add_argument(
        '--device', type=str, default=None,
                        help="Specify the GPUs for calling default: 'cuda:auto'")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
