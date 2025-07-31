import sys
import os
import copy
import csv
import tensorflow as tf
import logging
from time import time
import numpy as np

from argparse import ArgumentParser, SUPPRESS
from concurrent.futures import ProcessPoolExecutor, as_completed

from shared.utils import str2bool, log_error, log_warning, str_none, get_header
from shared.interval_tree import bed_tree_from
from clair3.CallVariants import OutputConfig, output_utilties_from, batch_output

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
logging.basicConfig(format='%(message)s', level=logging.INFO)


def Run(args):
    threads = args.threads
    os.environ["OMP_NUM_THREADS"] = f"{threads}"
    os.environ["OPENBLAS_NUM_THREADS"] = f"{threads}"
    os.environ["MKL_NUM_THREADS"] = f"{threads}"
    os.environ["NUMEXPR_NUM_THREADS"] = f"{threads}"

    # https://stackoverflow.com/a/54832345
    tf.config.threading.set_intra_op_parallelism_threads(threads)
    tf.config.threading.set_inter_op_parallelism_threads(2 if threads > 1 else 1)

    if not args.use_gpu:
        # https://www.intel.com/content/www/us/en/developer/articles/technical/maximize-tensorflow-performance-on-cpu-considerations-and-recommendations-for-inference.html
        os.environ["TF_ENABLE_ONEDNN_OPTS"] = "1"

    global test_pos
    test_pos = None
    global param
    if args.pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param

    if args.enable_long_indel:
        maximum_variant_length_that_need_infer = param.maximum_variant_length_that_need_infer_include_long_indel
    else:
        maximum_variant_length_that_need_infer = param.maximum_variant_length_that_need_infer

    output_config = OutputConfig(
        is_show_reference=args.showRef,
        is_debug=args.debug,
        is_haploid_precise_mode_enabled=args.haploid_precise,
        is_haploid_sensitive_mode_enabled=args.haploid_sensitive,
        is_output_for_ensemble=args.output_for_ensemble,
        quality_score_for_pass=args.qual,
        tensor_fn=args.tensor_fn,
        input_probabilities=args.input_probabilities,
        add_indel_length=args.add_indel_length,
        gvcf=args.gvcf,
        pileup=args.pileup,
        enable_long_indel=args.enable_long_indel,
        maximum_variant_length_that_need_infer=maximum_variant_length_that_need_infer,
        keep_iupac_bases=args.keep_iupac_bases
    )
    output_utilities = output_utilties_from(
        sample_name=args.sampleName,
        is_debug=args.debug,
        is_output_for_ensemble=args.output_for_ensemble,
        reference_file_path=args.ref_fn,
        output_file_path=args.call_fn,
        output_probabilities=args.output_probabilities,
        cmd_fn=args.cmd_fn
    )

    call_variants_from_cffi(args=args, output_config=output_config, output_utilities=output_utilities)

def tensor_generator_for_chunk(gen_cls, args, batch_size=50):
    if args.output_tensor_can_fn_list is not None:
        f_list = open(args.output_tensor_can_fn_list, 'r').read()
        f_list = f_list.strip().split('\n')
        import numpy as np
        for f in f_list:
            if f == '': continue
            parent_dir = os.path.dirname(args.output_tensor_can_fn_list)
            tensor = np.load(os.path.join(parent_dir, f +'.npy'))
            pos_alt = os.path.join(parent_dir, f + '.info')
            pos_alt_info = open(pos_alt, 'r').read().strip().split('\n')
            all_alt_info = []
            all_position = []
            for pos_alt in pos_alt_info:
                pos_alt = pos_alt.split('\t')
                all_position.append(pos_alt[0])
                all_alt_info.append(pos_alt[1])
            total_data = len(tensor)
            assert total_data == len(all_alt_info)
            assert total_data == len(all_position)
            batch_size = batch_size
            total_chunk = total_data // batch_size if total_data % batch_size == 0 else total_data // batch_size + 1
            for chunk_id in range(total_chunk):
                chunk_start = chunk_id * batch_size
                chunk_end = (chunk_id + 1) * batch_size if chunk_id < total_chunk - 1 else total_data
                yield (tensor[chunk_start:chunk_end],
                       all_position[chunk_start:chunk_end],
                       all_alt_info[chunk_start:chunk_end])
    else:
        tensor, all_position, all_alt_info = gen_cls(args)
        if tensor is None:
            return
        total_data = len(tensor)
        assert total_data == len(all_alt_info)
        assert total_data == len(all_position)
        batch_size = batch_size
        total_chunk = total_data // batch_size if total_data % batch_size == 0 else total_data // batch_size + 1
        for chunk_id in range(total_chunk):
            chunk_start = chunk_id * batch_size
            chunk_end = (chunk_id + 1) * batch_size if chunk_id < total_chunk - 1 else total_data
            yield (tensor[chunk_start:chunk_end],
                   all_position[chunk_start:chunk_end],
                   all_alt_info[chunk_start:chunk_end])

def tensor_generator_for_unchunked(gen_cls, args, batch_size=200):
    if args.bed_fn:
        tree = bed_tree_from(args.bed_fn)
        for contig in tree.values():
            contig.merge_overlaps()
        chunks = [(ctg_name, I.begin, I.end) for ctg_name, ctg_tree in tree.items() for I in ctg_tree.all_intervals]
    else:
        fai = f"{args.ref_fn}.fai"
        chunks = [(ctg, 0, ctg_len) for ctg, ctg_len, *_ in csv.reader(open(fai), delimiter='\t')]
    for ctg, start, end in chunks:
        temp_args = copy.copy(args)
        temp_args.ctgName = ctg
        temp_args.ctgStart = int(start)
        temp_args.ctgEnd = int(end)
        for batch in tensor_generator_for_chunk(gen_cls, temp_args, batch_size=batch_size):
            yield batch




def batch_output_worker(position, alt_info_list, shm_name, shape, dtype,
                        output_config, output_utilities):
    from multiprocessing import shared_memory
    try:
        shm = shared_memory.SharedMemory(name=shm_name)
        Y = np.ndarray(shape, dtype=dtype, buffer=shm.buf).copy()

        result = batch_output(position, alt_info_list, Y, output_config, output_utilities)

        return result
    finally:
        if 'shm' in locals():
            shm.close()
            shm.unlink()

def call_variants_from_cffi(args, output_config, output_utilities):
    if args.pileup:
        from preprocess.CreateTensorPileupFromCffi import CreateTensorPileup as CT
    else:
        from preprocess.CreateTensorFullAlignmentFromCffi import CreateTensorFullAlignment as CT

    if args.output_tensor_can_fn_list is None and args.use_gpu:
        CT(args)
        return

    use_triton_gpu = args.use_triton_gpu
    if use_triton_gpu:
        import tritonclient.grpc as tritongrpcclient
        server_url = 'localhost:8001'
        try:
            triton_client = tritongrpcclient.InferenceServerClient(
            url=server_url,
            verbose=False
        )
        except Exception as e:
            print("channel creation failed: " + str(e))
            sys.exit()
    elif args.use_gpu:
        gpus = tf.config.experimental.list_physical_devices('GPU')
        for gpu in gpus:
            tf.config.experimental.set_virtual_device_configuration(gpu, [
                tf.config.experimental.VirtualDeviceConfiguration(memory_limit=args.max_gpu_memory)])

        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id) if args.gpu_id else ""
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    global param
    if args.pileup:
        import shared.param_p as param
        if use_triton_gpu:
            model_name = 'pileup'
            input_dtype = 'INT32'
        else:
            from clair3.model import Clair3_P
            m = Clair3_P(add_indel_length=args.add_indel_length, predict=True)
    else:
        import shared.param_f as param
        if use_triton_gpu:
            model_name = 'alignment'
            input_dtype = 'INT8'
        else:
            from clair3.model import Clair3_F
            m = Clair3_F(add_indel_length=args.add_indel_length, predict=True)

    if not use_triton_gpu:
        m.load_weights(args.chkpnt_fn)
    args.output_file = open(args.call_fn, 'w') if args.call_fn != 'PIPE' else sys.stdout

    header_str = get_header(reference_file_path=args.ref_fn, cmd_fn=args.cmd_fn, sample_name=args.sampleName)
    header_str += '\n'
    args.output_file.write(header_str)

    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    full_alignment_mode = not args.pileup

    logging.info("Calling variants ...")
    variant_call_start_time = time()

    batch_output_method = batch_output
    total = 0

    if not args.use_gpu:
        batch_size = param.predictBatchSize * args.threads
    else:
        # increase bs in gpu for speed up
        batch_size = param.predictBatchSize * 5
    if args.ctgName == 'None' and args.pileup:
        batch_gen = tensor_generator_for_unchunked(CT, args, batch_size=batch_size)
    else:
        batch_gen = tensor_generator_for_chunk(CT, args, batch_size=batch_size)

    if not args.use_gpu:
        for (X, position, alt_info_list) in batch_gen:
            total += len(X)
            if args.pileup:
                for alt_idx, alt_info in enumerate(alt_info_list):
                    depth = int(alt_info.split('-', maxsplit=1)[0])
                    max_depth = param.max_depth_dict[args.platform]
                    # for extreme high coverage data, make sure we could have a truncated coverage
                    if depth > 0 and depth > max_depth * 1.5:
                        scale_factor = depth / max_depth
                        X[alt_idx] = X[alt_idx] / scale_factor

            if use_triton_gpu:
                inputs = [];
                outputs = []
                inputs.append(tritongrpcclient.InferInput('input_1', X.shape, input_dtype))
                outputs.append(tritongrpcclient.InferRequestedOutput('output_1'))
                inputs[0].set_data_from_numpy(X)
                results = triton_client.infer(model_name=model_name, inputs=inputs, outputs=outputs)
                Y = results.as_numpy('output_1')
            else:
                Y = m.predict_on_batch(X)

            output = batch_output_method(position, alt_info_list, Y, output_config, output_utilities, args)
    else:
        from multiprocessing import shared_memory
        with ProcessPoolExecutor(max_workers=args.cpu_threads) as executor:
            futures = []
            for X, position, alt_info_list in batch_gen:
                total += len(X)
                Y = m.predict_on_batch(X)
                shm = shared_memory.SharedMemory(create=True, size=Y.nbytes)
                shared_Y = np.ndarray(Y.shape, dtype=Y.dtype, buffer=shm.buf)
                np.copyto(shared_Y, Y)

                future = executor.submit(
                    batch_output_worker,
                    position, alt_info_list,
                    shm.name, Y.shape, Y.dtype,
                    output_config, output_utilities
                )
                futures.append(future)

                shm.close()

            for future in as_completed(futures):
                try:
                    result = future.result()
                    args.output_file.write(result)
                except Exception as e:
                    print(f"[Error] Task failed: {e}", flush=True)

    if chunk_id is not None:
        logging.info("Total processed positions in {} (chunk {}/{}) : {}".format(args.ctgName, chunk_id+1, chunk_num, total))
    elif full_alignment_mode:
        try:
            chunk_infos = args.call_fn.split('.')[-2]
            c_id, c_num = chunk_infos.split('_')
            c_id = int(c_id) + 1 # 0-index to 1-index
            logging.info("Total processed positions in {} (chunk {}/{}) : {}".format(args.ctgName, c_id, c_num, total))
        except:
            logging.info("Total processed positions {}: {}".format(args.ctgName if args.ctgName else "", total))
    else:
        logging.info("Total processed positions {}: {}".format(args.ctgName if args.ctgName else "", total))

    if full_alignment_mode and total == 0:
        logging.info(log_warning("[WARNING] No full-alignment output for file {}/{}".format(args.ctgName, args.call_fn)))

    logging.info("Total time elapsed: %.2f s" % (time() - variant_call_start_time))

    args.output_file.close()
    # remove file if on variant in output
    if os.path.exists(args.call_fn):
        for row in open(args.call_fn, 'r'):
            if row[0] != '#':
                return
        logging.info("[INFO] No vcf output for file {}, remove empty file".format(args.call_fn))
        os.remove(args.call_fn)


def main():
    parser = ArgumentParser(description="Call variants using a trained model and tensors of candidate variants")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
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

    # options for debug purpose
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help="DEBUG: Use GPU for calling. Speed up is mostly insignficiant. Only use this for building your own pipeline")

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

    ## Use bin file from pytables to speed up calling.
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

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--gpu_id', type=str, default=None,
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--max_gpu_memory', type=int, default=5000,
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--use_triton_gpu', type=str2bool, default=False,
                        help="DEBUG: Use GPU for calling. Speed up is mostly insignficiant. Only use this for building your own pipeline")

    parser.add_argument('--output_tensor_can_fn_list', type=str, default=None,
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--cpu_threads', type=int, default=4,
                        help="How many threads for python to use for parallelization default: %(default)f")

    parser.add_argument('--output_file', type=str, default=None,
                        help="How many threads for python to use for parallelization default: %(default)f")


    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()