import os
import shlex
import logging
import numpy as np
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import libclair3
import shared.param_f as param
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, str2bool, vcf_candidates_from
from shared.interval_tree import bed_tree_from

logging.basicConfig(format='%(message)s', level=logging.INFO)
no_of_positions = param.no_of_positions
flanking_base_num = param.flankingBaseNum
channel_size = param.channel_size


def CreateTensorFullAlignment(args):

    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    full_aln_regions = args.full_aln_regions
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    bam_file_path = args.bam_fn
    extend_bp = param.extend_bp
    platform = args.platform
    phased_vcf_fn = args.phased_vcf_fn
    min_mapping_quality = args.minMQ
    min_base_quality = args.minBQ
    enable_long_indel = args.enable_long_indel
    extend_bed = file_path_from(args.extend_bed)
    is_extend_bed_file_given = extend_bed is not None
    confident_bed_fn = file_path_from(args.bed_fn)
    is_confident_bed_file_given = confident_bed_fn is not None

    # we would't haplotag reads if --no_phasing_for_fa option is enabled
    need_haplotagging = args.no_phasing_for_fa is not True
    candidates_set = set()
    matrix_depth = param.matrix_depth_dict[platform]
    max_indel_length = param.maximum_variant_length_that_need_infer if not enable_long_indel else param.maximum_variant_length_that_need_infer_include_long_indel

    if full_aln_regions:

        """
        If given full alignment bed regions, all candidate positions will be directly selected from each row, define as 
        'ctg start end', where 0-based center position is the candidate for full alignment calling.
        if 'need_haplotagging' option enables, full alignment bed regions will also include nearby heterozygous snp candidates for reads
        haplotag, which is faster than whatshap haplotag with more memory occupation.
        """

        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (full_aln_regions)))
        candidate_file_path_output = candidate_file_path_process.stdout

        ctg_start, ctg_end = float('inf'), 0
        for row in candidate_file_path_output:
            row = row.rstrip().split('\t')
            if row[0] != ctg_name: continue
            position = int(row[1]) + 1
            end = int(row[2]) + 1
            ctg_start = min(position, ctg_start)
            ctg_end = max(end, ctg_end)

            if len(row) > 3:  # hete snp positions
                center_pos = position + extend_bp + 1
                ref_base, alt_base, genotype, phase_set = row[3].split('-')
            else:
                center = position + (end - position) // 2 - 1
                candidates_set.add(center)

        candidate_file_path_output.close()
        candidate_file_path_process.wait()

    variant_list = []
    if need_haplotagging and phased_vcf_fn and os.path.exists(phased_vcf_fn):
        # if need_haplotagging option enables, scan the phased vcf file and store the heterozygous SNP candidates from each phase set
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (phased_vcf_fn)))
        for row in unzip_process.stdout:
            row = row.rstrip()
            if row[0] == '#':
                continue
            columns = row.strip().split('\t')
            contig_name = columns[0]
            if ctg_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            genotype_info = columns[9].split(':')
            genotype, phase_set = genotype_info[0], genotype_info[-1]
            if '|' not in genotype:  # unphasable
                continue
            genotype = ('1' if genotype == '0|1' else '2')

            # use a C Variant struct to store all phased infos
            variant_list.append(libclair3.ffi.new("struct Variant *", [pos-1, ref_base.encode(), alt_base.encode(), int(genotype), int(phase_set)]))

        variant_num = len(variant_list)
        Variants = libclair3.ffi.new("struct Variant *[]", variant_list)
    else:
        Variants = libclair3.ffi.new("struct Variant *[]", 1)
        variant_num = 0

    # 1-index to 0-index
    candidates_list = sorted(list(set([item-1 for item in candidates_set if item >= ctg_start and item <= ctg_end])))

    region_str = '{}:{}-{}'.format(ctg_name, ctg_start, ctg_end).encode()
    candidate_num = len(candidates_list)

    candidates = libclair3.ffi.new("size_t [{}]".format(candidate_num), candidates_list)

    fa_data = libclair3.lib.calculate_clair3_full_alignment(region_str, bam_file_path.encode(), fasta_file_path.encode(),
                                                      Variants, variant_num, candidates, candidate_num, need_haplotagging,
                                                            min_mapping_quality, min_base_quality, matrix_depth, max_indel_length)

    # use np buffer to get the matrix
    matrix_depth = param.matrix_depth_dict[platform]
    ffi = libclair3.ffi
    _dtype = np.int8
    size_sizet = np.dtype(_dtype).itemsize
    np_fa_data = np.frombuffer(ffi.buffer(
        fa_data.matrix, size_sizet * matrix_depth * no_of_positions * channel_size * candidate_num),
        dtype=_dtype
    ).reshape(candidate_num, matrix_depth, no_of_positions, channel_size).copy()


    all_position_info, all_alt_info = [], []
    for idx in range(candidate_num):
        # decode the C char* to python string
        alt_info_string = ffi.string(fa_data.all_alt_info[idx]).decode('utf8', 'ignore')
        alt_info = alt_info_string.rstrip().split('-')
        pos, depth, center_ref_base, alt = alt_info[:4]
        all_position_info.append(ctg_name + ':' + pos + ':' + center_ref_base)
        all_alt_info.append(depth + '-' + alt)

    libclair3.lib.destroy_fa_data(fa_data)

    return np_fa_data, all_position_info, all_alt_info


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

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
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

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
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # options for advanced users
    parser.add_argument('--minCoverage', type=int, default=param.min_coverage,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$minMQ are filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$minBQ are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=param.max_depth,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--phasing_window_size', type=int, default=param.phasing_window_size,
                        help="DEBUG: The window size for read phasing")

    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--indel_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=5, type=int,
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

    ## Use heterozygous SNP variants in phased vcf file for haplotaging
    parser.add_argument('--phased_vcf_fn', type=str, default=None,
                        help=SUPPRESS)
    ## Apply no phased data in training. Only works in data training, default: False
    parser.add_argument('--add_no_phasing_data_training', action='store_true',
                        help=SUPPRESS)

    ## Output representation unification infos, which refines training labels
    parser.add_argument('--unify_repre', action='store_true',
                        help=SUPPRESS)

    ## Path of representation unification output
    parser.add_argument('--unify_repre_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Provide the regions to be included in full-alignment based calling
    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help=SUPPRESS)

    ## Use Clair3's own phasing module for read level phasing when creating tensor, compared to using Whatshap, speed is faster but has higher memory footprint, default: False
    parser.add_argument('--need_haplotagging', action='store_true',
                        help=SUPPRESS)

    ## Apply read realignment for illumina platform. Greatly boost indel performance in trade of running time
    parser.add_argument('--need_realignment', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    CreateTensorFullAlignment(args)


if __name__ == "__main__":
    main()