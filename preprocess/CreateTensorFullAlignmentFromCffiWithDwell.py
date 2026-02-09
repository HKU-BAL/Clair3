"""
CreateTensorFullAlignmentFromCffi with Dwell Time Support

This is an enhanced version of CreateTensorFullAlignmentFromCffi.py that includes
support for nanopore dwell time extraction from mv tags.

Key differences:
1. Added enable_dwell_time parameter
2. Adjusts channel_size based on dwell time option (8 or 9 channels)
3. Passes enable_dwell_time to C function

Usage:
    python CreateTensorFullAlignmentFromCffiWithDwell.py \
        --bam_fn input.bam \
        --ref_fn reference.fa \
        --tensor_can_fn output.npy \
        --full_aln_regions regions.bed \
        --ctgName chr1 \
        --enable_dwell_time
"""

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
    enable_dwell_time = args.enable_dwell_time
    extend_bed = file_path_from(args.extend_bed)
    is_extend_bed_file_given = extend_bed is not None
    confident_bed_fn = file_path_from(args.bed_fn)
    is_confident_bed_file_given = confident_bed_fn is not None

    # Adjust channel size based on dwell time option
    channel_size = param.channel_size + 1 if enable_dwell_time else param.channel_size

    # we would't haplotag reads if --no_phasing_for_fa option is enabled
    need_haplotagging = args.no_phasing_for_fa is not True
    candidates_set = set()
    matrix_depth = param.matrix_depth_dict[platform]
    max_indel_length = param.maximum_variant_length_that_need_infer if not enable_long_indel else param.maximum_variant_length_that_need_infer_include_long_indel

    all_position_info, all_alt_info = [], []

    if full_aln_regions:

        """
        If given full alignment bed regions, all candidate positions will be directly selected from each row, define as 
        'ctg start end', where 0-based center position is the candidate for full alignment calling.
        if 'need_haplotagging' option enables, full alignment bed regions will also include nearby heterozygous snp candidates for reads
        haplotag, which is faster than whatshap haplotag with more memory occupation.
        """

        ctg_start, ctg_end = float('inf'), 0
        with open(full_aln_regions, 'r') as f:
            for row in f:
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
                    if position == 1:
                        center = end - flanking_base_num - 2
                    else:
                        center = position + (end - position) // 2 - 1
                    candidates_set.add(center)

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

    if ctg_start is None or ctg_end is None:
        return [], all_position_info, all_alt_info
    # 1-index to 0-index
    candidates_list = sorted(list(set([item-1 for item in candidates_set if item >= ctg_start and item <= ctg_end])))

    region_str = '{}:{}-{}'.format(ctg_name, ctg_start, ctg_end).encode()
    candidate_num = len(candidates_list)

    candidates = libclair3.ffi.new("size_t [{}]".format(candidate_num), candidates_list)

    # Call C function with enable_dwell_time parameter
    fa_data = libclair3.lib.calculate_clair3_full_alignment(
        region_str, 
        bam_file_path.encode(), 
        fasta_file_path.encode(),
        Variants, 
        variant_num, 
        candidates, 
        candidate_num, 
        need_haplotagging,
        min_mapping_quality, 
        min_base_quality, 
        matrix_depth, 
        max_indel_length,
        enable_dwell_time  # NEW: Enable dwell time extraction
    )

    # use np buffer to get the matrix with correct channel size
    matrix_depth = param.matrix_depth_dict[platform]
    ffi = libclair3.ffi
    _dtype = np.int8
    size_sizet = np.dtype(_dtype).itemsize
    np_fa_data = np.frombuffer(ffi.buffer(
        fa_data.matrix, size_sizet * matrix_depth * no_of_positions * channel_size * candidate_num),
        dtype=_dtype
    ).reshape(candidate_num, matrix_depth, no_of_positions, channel_size).copy()

    for idx in range(candidate_num):
        # decode the C char* to python string
        alt_info_string = ffi.string(fa_data.all_alt_info[idx]).decode('utf8', 'ignore')
        alt_info = alt_info_string.rstrip().split('-')
        pos, depth, center_ref_base, alt = alt_info[:4]
        all_position_info.append(ctg_name + ':' + pos + ':' + center_ref_base)
        all_alt_info.append(depth + '-' + alt)

    libclair3.lib.destroy_fa_data(fa_data)
    
    if args.tensor_can_fn != "PIPE":
        if len(all_position_info) == 0:
            return None, None, None
        np.save(args.tensor_can_fn, np.array(np_fa_data, dtype=np.int8))
        with open(args.tensor_can_fn + '.info', 'w') as info_file:
            for idx, pos_info in enumerate(all_position_info):
                info_file.write(f"{pos_info}\t{all_alt_info[idx]}\n")
        print("[INFO] Total processed positions in {} : {}".format(args.ctgName, len(all_position_info)))
        print(f"[INFO] Output matrix shape: {np_fa_data.shape} (with {channel_size} channels)")
        if enable_dwell_time:
            print("[INFO] Dwell time channel (channel 8) is enabled")
        return None, None, None

    return np_fa_data, all_position_info, all_alt_info


def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment with dwell time support")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam", required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions")

    parser.add_argument('--extend_bed', type=str, default=None,
                        help="Extend the regions in the --bed_fn")

    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help="Provide the regions to be included in full-alignment based calling")

    parser.add_argument('--phased_vcf_fn', type=str, default=None,
                        help="Phased VCF file for haplotagging")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="Minimum mapping quality, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="Minimum base quality, default: %(default)d")

    parser.add_argument('--no_phasing_for_fa', action='store_true',
                        help="Disable phasing for full alignment")

    parser.add_argument('--enable_long_indel', action='store_true',
                        help="Enable long indel calling")

    parser.add_argument('--enable_dwell_time', action='store_true',
                        help="Enable dwell time extraction from mv tag (adds channel 8)")

    args = parser.parse_args()

    CreateTensorFullAlignment(args)


if __name__ == "__main__":
    main()

