import sys
import logging
import queue
import concurrent.futures
import numpy as np

from argparse import ArgumentParser, SUPPRESS
from contextlib import contextmanager

import libclair3
import shared.param_p as param
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import file_path_from, IUPAC_base_to_num_dict as BASE2NUM, str2bool, vcf_candidates_from
from preprocess.medaka_utils import Region

logging.getLogger().setLevel(logging.INFO)

flanking_base_num = param.flankingBaseNum
no_of_positions = 2 * flanking_base_num + 1
channel = param.channel
channel_size = len(channel)


def pileup_counts_clair3(
        region, bam, fasta, min_depth, min_snp_af, min_indel_af, min_mq, call_snp_only, max_indel_length, \
        max_depth, gvcf=False, region_split=100000, workers=1):
    """Create pileup counts feature array for region.

    :param region: `medaka.common.Region` object
    :param bam: .bam file with alignments.
    :param dtype_prefixes: prefixes for query names which to separate counts.
        If `None` (or of length 1), counts are not split.
    :param region_split: largest region to process in single thread.
        Regions are processed in parallel and stitched before being returned.
    :param workers: worker threads for calculating pileup.
    :param tag_name: two letter tag name by which to filter reads.
    :param tag_value: integer value of tag for reads to keep.
    :param keep_missing: whether to keep reads when tag is missing.
    :param num_qstrat: number of layers for qscore stratification.
    :param weibull_summation: use a Weibull partial-counts approach,
        requires 'WL' and 'WK' float-array tags.

    :returns: iterator of tuples
        (pileup counts array, reference positions, insertion positions)
        Multiple chunks are returned if there are discontinuities in
        positions caused e.g. by gaps in coverage.
    """
    lib = libclair3.lib
    featlenclair3 = lib.featlenclair3
    bam = BAMHandler(bam)

    def _process_region(reg):
        # ctg start is 1-based, medaka.common.Region object is 0-based
        region_str = '{}:{}-{}'.format(reg.ref_name, max(0, reg.start-1), reg.end)
        if isinstance(bam, BAMHandler):
            bam_handle = bam
        else:
            bam_handle = BAMHandler(bam)
        with bam_handle.borrow() as fh:
            counts = lib.calculate_clair3_pileup(
                region_str.encode(), fh, fasta.encode(), min_depth, min_snp_af, min_indel_af, min_mq, max_indel_length, call_snp_only, max_depth, gvcf)

        np_counts, positions, alt_info_string_list, gvcf_output = _plp_data_to_numpy(
            counts, featlenclair3, gvcf=gvcf)

        alt_info_list = []
        for alt_info in alt_info_string_list:
            alt_info = alt_info.split('-')
            # skip mainly because candidate length is larger than maximum indel length
            if len(alt_info) < 4:
                continue
            pos, depth, center_ref_base, alt = alt_info[:4]
            alt_info_list.append((int(pos), reg.ref_name + ':' + pos + ':' + center_ref_base, depth + '-' + alt))

        lib.destroy_plp_data(counts, gvcf)
        return np_counts, positions, alt_info_list, gvcf_output

    # we found that split into small chunk would lead to some missing truths,
    # the candidates cross two negbouring small chunks
    region_split = region.end - region.start
    regions = region.split(region_split, fixed_size=False)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        results = executor.map(_process_region, regions)
        chunk_results, all_alt_info_list, gvcf_output = __enforce_pileup_chunk_contiguity(results)
    return chunk_results, all_alt_info_list, gvcf_output


class BAMHandler(object):
    """Opening of BAM file handles and indices."""

    def __init__(self, bam, size=16):
        """Initialise a pool of HTSlib filehandles."""
        # note: the default size here is set to match the default
        #       `bam_workers` of prediction.DataLoader and `workers`
        #       of features.pileup_counts, such that this class
        #       should never block computations
        self.bam = bam
        self._pool = queue.Queue(size)

        lib, ffi = libclair3.lib, libclair3.ffi
        for _ in range(size):
            fset = ffi.gc(
                lib.create_bam_fset(self.bam.encode()),
                self._destroy_fset)
            self._pool.put(fset)

    @contextmanager
    def borrow(self):
        """Borrow a BAM file handle and index set."""
        fset = self._pool.get()
        try:
            yield fset
        finally:
            self._pool.put(fset)

    def encode(self):
        """Return bare path encoded to bytes.

        For legacy compatibility only.
        """
        return self.bam.encode()

    def _destroy_fset(self, fset):
        libclair3.lib.destroy_bam_fset(fset)


def _plp_data_to_numpy(plp_data, n_rows, gvcf=False):
    """Create numpy representation of feature data.

    Copy the feature matrix and alignment column names from a
    `plp_data` structure returned from C library function calls.

    :param plp_data: a cffi proxy to a `plp_data*` pointer
    :param nrows: the number of rows in the plp_data.matrix (the number
        of elements in the feature per pileup column).

    :returns: pileup counts numpy array, reference positions

    """
    ffi = libclair3.ffi
    size_sizet = np.dtype(int).itemsize
    _dtype = int
    np_counts = np.frombuffer(ffi.buffer(
        plp_data.matrix, size_sizet * plp_data.n_cols * n_rows),
        dtype=_dtype
    ).reshape(plp_data.n_cols, n_rows).copy()

    alt_info_string_list = []
    gvcf_output = []
    candidates_num = plp_data.candidates_num
    # decode all alternative information, position-depth-reference_base-alt_info
    for i in range(candidates_num):
        alt_info_string = ffi.string(plp_data.all_alt_info[i]).decode('utf8', 'ignore').rstrip()
        alt_info_string_list.append(alt_info_string)

    if gvcf:
        gvcf_pos_ref_count = np.frombuffer(ffi.buffer(
            plp_data.pos_ref_count, size_sizet * plp_data.buffer_cols),
            dtype=_dtype
        ).reshape(plp_data.buffer_cols).copy()
        gvcf_pos_total_count = np.frombuffer(ffi.buffer(
            plp_data.pos_total_count, size_sizet * plp_data.buffer_cols),
            dtype=_dtype
        ).reshape(plp_data.buffer_cols).copy()
        gvcf_output = [gvcf_pos_ref_count, gvcf_pos_total_count]

    positions = np.empty(plp_data.n_cols, dtype=[
        ('major', int), ('minor', int)])
    np.copyto(
        positions['major'], np.frombuffer(
            ffi.buffer(plp_data.major, size_sizet * plp_data.n_cols),
            dtype=_dtype))
    np.copyto(
        positions['minor'],
        np.frombuffer(ffi.buffer(
            plp_data.minor, size_sizet * plp_data.n_cols), dtype=_dtype))
    return np_counts, positions, alt_info_string_list, gvcf_output


def __enforce_pileup_chunk_contiguity(pileups):
    """Split and join ordered pileup chunks to ensure contiguity.

    :param pileups: iterable of (counts, pileups) as constructed by
        `_plp_data_to_numpy`.

    :returns: a list of reconstituted (counts, pileups) where discontinuities
        in the inputs cause breaks and abutting inputs are joined.

    """
    split_results = list()
    all_alt_info_list = list()
    # First pass: need to check for discontinuities within chunks,
    # these show up as >1 changes in the major coordinate
    for counts, positions, alt_info_list, gvcf_output in pileups:
        move = np.ediff1d(positions['major'])
        gaps = np.where(move > 1)[0] + 1
        all_alt_info_list += alt_info_list
        if len(gaps) == 0:
            split_results.append((counts, positions))
        else:
            start = 0
            for i in gaps:
                split_results.append((counts[start:i], positions[start:i]))
                start = i
            split_results.append((counts[start:], positions[start:]))

    # Second pass: stitch abutting chunks together, anything not neighbouring
    # is kept separate whether it came from the same chunk originally or not
    def _finalize_chunk(c_buf, p_buf):
        chunk_counts = np.concatenate(c_buf)
        chunk_positions = np.concatenate(p_buf)
        return chunk_counts, chunk_positions

    counts_buffer, positions_buffer = list(), list()
    chunk_results = list()
    last = None
    for counts, positions in split_results:
        if len(positions) == 0:
            continue
        first = positions['major'][0]
        # should be last -first == 1?
        if len(counts_buffer) == 0 or first - last == 1:
            # new or contiguous
            counts_buffer.append(counts)
            positions_buffer.append(positions)
            last = positions['major'][-1]
        else:
            # discontinuity
            chunk_results.append(_finalize_chunk(
                counts_buffer, positions_buffer))
            counts_buffer = [counts]
            positions_buffer = [positions]
            last = positions['major'][-1]
    if len(counts_buffer) != 0:
        chunk_results.append(_finalize_chunk(counts_buffer, positions_buffer))
    return chunk_results, all_alt_info_list, gvcf_output


def CreateTensorPileup(args):
    """
    Create pileup tensor for pileup model training or calling.
    Use slide window to scan the whole candidate regions, keep all candidates over specific minimum allelic frequency
    and minimum depth, use samtools mpileup to store pileup info for pileup tensor generation. Only scan candidate
    regions once, we could directly get all variant candidates directly.
    """
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    bam_file_path = args.bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    minimum_snp_af_for_candidate = args.snp_min_af
    minimum_indel_af_for_candidate = args.indel_min_af
    min_coverage = args.minCoverage
    min_mapping_quality = args.minMQ
    platform = args.platform

    vcf_fn = file_path_from(args.vcf_fn)
    is_known_vcf_file_provided = vcf_fn is not None
    confident_bed_fn = file_path_from(args.extend_bed)
    is_confident_bed_file_given = confident_bed_fn is not None
    extend_bed = file_path_from(args.extend_bed)
    is_extend_bed_file_given = extend_bed is not None
    fast_mode = args.fast_mode
    call_snp_only = args.call_snp_only
    enable_long_indel = args.enable_long_indel
    # 1-based regions [start, end] (start and end inclusive)
    tree, bed_start, bed_end = bed_tree_from(bed_file_path=extend_bed,
                                             contig_name=ctg_name,
                                             return_bed_region=True)

    fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')

    fast_mode = platform == 'ont' and fast_mode
    minimum_snp_af_for_candidate = max(minimum_snp_af_for_candidate, param.min_af_dict[platform]) if fast_mode else minimum_snp_af_for_candidate
    min_coverage = max(min_coverage, 4) if fast_mode else min_coverage
    max_indel_length = param.maximum_variant_length_that_need_infer if not enable_long_indel else param.maximum_variant_length_that_need_infer_include_long_indel

    if not is_confident_bed_file_given and chunk_id is not None:
        contig_length = 0
        with open(fai_fn, 'r') as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")

                contig_name = columns[0]
                if contig_name != ctg_name:
                    continue
                contig_length = int(columns[1])
        chunk_size = contig_length // chunk_num + 1 if contig_length % chunk_num else contig_length // chunk_num
        ctg_start = chunk_size * chunk_id  # 0-base to 1-base
        ctg_end = ctg_start + chunk_size

    if is_confident_bed_file_given and chunk_id is not None:
        chunk_size = (bed_end - bed_start) // chunk_num + 1 if (bed_end - bed_start) % chunk_num else (bed_end - bed_start) // chunk_num
        ctg_start = bed_start + 1 + chunk_size * chunk_id  # 0-base to 1-base
        ctg_end = ctg_start + chunk_size

    if is_known_vcf_file_provided and chunk_id is not None:
        known_variants_list = vcf_candidates_from(vcf_fn=vcf_fn, contig_name=ctg_name)
        total_variants_size = len(known_variants_list)
        chunk_variants_size = total_variants_size // chunk_num if total_variants_size % chunk_num == 0 else total_variants_size // chunk_num + 1
        chunk_start_pos = chunk_id * chunk_variants_size
        known_variants_set = set(known_variants_list[chunk_start_pos: chunk_start_pos + chunk_variants_size])
        if len(known_variants_set) == 0:
            return [], [], []
        ctg_start, ctg_end = min(known_variants_set), max(known_variants_set)

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    if is_ctg_range_given:
        ctg_start = max(1, ctg_start)
        extend_start = max(1, ctg_start - no_of_positions)
        extend_end = ctg_end + no_of_positions

    region_str = "{}:{}-{}".format(ctg_name, extend_start, extend_end)
    region = Region.from_string(region_str)

    confident_bed_tree = bed_tree_from(bed_file_path=confident_bed_fn, contig_name=ctg_name, bed_ctg_start=extend_start,
                                       bed_ctg_end=extend_end)

    if args.gvcf:
        from preprocess.utils import variantInfoCalculator
        nonVariantCaller = variantInfoCalculator(gvcfWritePath=args.temp_file_dir, ref_path=args.ref_fn,
                                                 bp_resolution=args.bp_resolution, ctgName=ctg_name,sample_name='.'.join(
                [args.sampleName, ctg_name, str(ctg_start), str(ctg_end)]), p_err=args.base_err,
                                                 gq_bin_size=args.gq_bin_size)

    chunk_result, all_alt_info_list, gvcf_output = pileup_counts_clair3(region,
                                                           bam=bam_file_path,
                                                           fasta=fasta_file_path,
                                                           min_depth=min_coverage,
                                                           min_snp_af=minimum_snp_af_for_candidate,
                                                           min_indel_af=minimum_indel_af_for_candidate,
                                                           min_mq=min_mapping_quality,
                                                           max_indel_length=max_indel_length,
                                                           call_snp_only=call_snp_only,
                                                           max_depth=param.max_depth,
                                                           gvcf=args.gvcf)

    # slice all candidates tensor according to the alternative information
    np_pileup_data, all_position_info, all_alt_info = [], [], []
    for idx, (pos, pos_info, alt_info) in enumerate(all_alt_info_list):
        pos = int(pos)
        pass_confident_bed = not is_confident_bed_file_given or is_region_in(tree=confident_bed_tree,
                                                                             contig_name=ctg_name,
                                                                             region_start=pos - 1,
                                                                             region_end=pos + 1)

        pass_vcf_region = not is_known_vcf_file_provided or (is_known_vcf_file_provided and pos in known_variants_set)

        if not pass_confident_bed or not pass_vcf_region:
            continue
        start, end = pos - flanking_base_num, pos + flanking_base_num + 1
        for result in chunk_result:
            if start - 1 >= result[1][0][0] and end <= result[1][-1][0]:
                offset = start - result[1][0][0] - 1
                tensor = result[0][offset: offset+no_of_positions]
                # mainly because no coverage in flanking windows
                if tensor.shape != (no_of_positions, channel_size):
                    continue
                # check any empty columns in flanking position, those columns with all zeros
                if np.sum(np.sum(tensor == 0, axis=1) == channel_size) > 0:
                    continue
                np_pileup_data.append(tensor)
                all_position_info.append(pos_info)
                all_alt_info.append(alt_info)
    np_pileup_data = np.array(np_pileup_data, dtype=np.int32)


    if args.gvcf:

        from shared.utils import reference_sequence_from, region_from
        samtools_execute_command = args.samtools
        ref_regions = []
        reference_start, reference_end = ctg_start - param.expandReferenceRegion, ctg_end + param.expandReferenceRegion
        reference_start = 1 if reference_start < 1 else reference_start
        ref_regions.append(region_from(ctg_name=ctg_name, ctg_start=reference_start, ctg_end=reference_end))
        reference_sequence = reference_sequence_from(
            samtools_execute_command=samtools_execute_command,
            fasta_file_path=fasta_file_path,
            regions=ref_regions
        )

        offset = 0 if ctg_start == 1 else 1
        empty_pileup_flag = False
        start = ctg_start - extend_start + offset
        end = ctg_end + 1 - extend_start + offset
        if sum(gvcf_output[1][start:end]) == 0:
            empty_pileup_flag = True
        for pos in range(ctg_start, ctg_end):
            if empty_pileup_flag:
                break
            ref_count = gvcf_output[0][pos - extend_start + offset]
            total_count = gvcf_output[1][pos - extend_start + offset]
            if pos -reference_start >= len(reference_sequence):
                continue
            reference_base = reference_sequence[pos-reference_start]
            if (ref_count == 0 and total_count == 0):
                cur_site_info = {'chr': ctg_name, 'pos': pos, 'ref': reference_base, 'n_total': 0, 'n_ref': 0}
                nonVariantCaller.make_gvcf_online(cur_site_info)
                continue

            cur_site_info = {'chr': ctg_name, 'pos': pos, 'ref': reference_base, 'n_total': total_count,
                             'n_ref': ref_count}
            nonVariantCaller.make_gvcf_online(cur_site_info)
        if len(nonVariantCaller.current_block) != 0:
            nonVariantCaller.write_to_gvcf_batch(nonVariantCaller.current_block, nonVariantCaller.cur_min_DP,
                                                 nonVariantCaller.cur_raw_gq)

        if empty_pileup_flag:
            nonVariantCaller.write_empty_pileup(ctg_name, ctg_start, ctg_end)
        nonVariantCaller.close_vcf_writer()

    return np_pileup_data, all_position_info, all_alt_info


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using pileup")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam", required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--min_af', type=float, default=0.08,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--snp_min_af', type=float, default=0.08,
                        help="Minimum snp allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--indel_min_af', type=float, default=0.15,
                        help="Minimum indel allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, default: %(default)s")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # options for advanced users
    parser.add_argument('--fast_mode', type=str2bool, default=False,
                        help="EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: %(default)s")

    parser.add_argument('--minCoverage', type=int, default=2,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$minMQ are filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$minBQ are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=param.max_depth,
                        help="EXPERIMENTAL: Maximum pileup depth to be processed. default: %(default)s")

    parser.add_argument('--call_snp_only', type=str2bool, default=False,
                        help="EXPERIMENTAL: Call candidates pass snp minimum AF only, ignore Indel candidates")

    # options for debug purpose
    parser.add_argument('--extend_bed', type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--temp_file_dir', type=str, default="./",
                        help="EXPERIMENTAL: The cache directory for storing temporary non-variant information if --gvcf is enabled, default: %(default)s")

    parser.add_argument('--indel_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    parser.add_argument('--base_err', default=param.base_err, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=param.gq_bin_size, type=int,
                        help='DEBUG: Default gq bin size for merge non-variant block in gvcf option, default: %(default)d')

    parser.add_argument('--bp_resolution', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")

    # options for internal process control
    ## Path to the 'zstd' compression
    parser.add_argument('--zstd', type=str, default=param.zstd,
                        help=SUPPRESS)

    ## Test in specific candidate position. Only for testing
    parser.add_argument('--test_pos', type=int, default=0,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    CreateTensorPileup(args)


if __name__ == "__main__":
    main()
