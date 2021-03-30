import sys
import shlex
import logging
from subprocess import PIPE
from os.path import isfile
from argparse import ArgumentParser
from collections import Counter, defaultdict

import shared.param_p as param
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool
from clair3.MakeGvcf import variantInfoCalculator

logging.getLogger().setLevel(logging.INFO)
BASES = set(list(BASE2NUM.keys()) + ["-"])
flanking_base_num = param.flankingBaseNum
sliding_window_size = no_of_positions = 2 * flanking_base_num + 1

BASE2NUMBER = dict(zip(
    "ACGTURYSWKMBDHVN-",
    (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 4)
))
channel = param.channel
channel_size = len(channel)
BASE2INDEX = dict(zip(channel, tuple(range(channel_size))))


def phredscore2raw_score(qual):
    return ord(qual) - 33


def evc_base_from(base):
    if base == 'N':
        return 'A'
    elif base == 'n':
        return 'a'
    elif base in 'ACGTacgt':
        return base
    elif base.isupper():
        return 'A'
    else:
        return 'a'


class CandidateStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def generate_tensor(pos, pileup_bases, reference_sequence, reference_start, reference_base, minimum_af_for_candidate,
                    minimum_snp_af_for_candidate, minimum_indel_af_for_candidate, platform, fast_mode):
    """
    Generate pileup input tensor
    pos: center position for pileup generation, default no_of_positions = flankingBaseNum + 1 + flankingBaseNum
    pileup_bases: pileup bases list of each read in specific candidate position from samtools mpileup 1.10
    reference_sequence: the whole reference sequence index by contig:start-end. 0-based.
    reference_base: upper reference base for cigar calculation.
    reference_start: 0 based reference start position for region querying.
    minimum_af_for_candidate: default minimum alleic frequency for candidate filtering, filter if below specific thredshold.
    """

    reference_base = evc_base_from(reference_base)
    pileup_tensor = [0] * channel_size
    base_idx = 0
    base_list = []
    alt_dict = defaultdict(int)
    pileup_dict = defaultdict(int)
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list.append(base + pileup_bases[base_idx: base_idx + advance])
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append(base)
        elif base == '^':  # start of a read, next character is mapping quality
            base_idx += 1
        # elif base == '$': # end of read with '$' symbol
        base_idx += 1
    base_counter = Counter(base_list)
    depth, max_ins_0, max_del_0, max_ins_1, max_del_1 = 0, 0, 0, 0, 0
    max_del_length = 0
    for key, count in base_counter.items():
        if key[0] == '+':
            alt_dict['I' + reference_base + key[1:].upper()] += count
            pileup_dict['I'] += count
            # two strand
            if key[1] in 'ACGTN*':
                pileup_tensor[BASE2INDEX["I"]] += count
                max_ins_0 = max(max_ins_0, count)
            else:
                pileup_tensor[BASE2INDEX["i"]] += count
                max_ins_1 = max(max_ins_1, count)
        elif key[0] == '-':
            del_base = reference_sequence[pos - reference_start + 1: pos - reference_start + len(key[1:]) + 1]
            alt_dict['D' + del_base] += count
            pileup_dict['D'] += count
            max_del_length = max(max_del_length, len(del_base))
            # two strand
            if key[1] in 'N*ACGT':
                pileup_tensor[BASE2INDEX["D"]] += count
                max_del_0 = max(max_del_0, count)
            else:
                pileup_tensor[BASE2INDEX["d"]] += count
        else:
            if key.upper() in 'ACGT':
                pileup_dict[key.upper()] += count
                depth += count
                if key.upper() != reference_base:
                    alt_dict['X' + key.upper()] += count
                pileup_tensor[BASE2INDEX[key]] += count
            elif key in '#*':
                pileup_tensor[BASE2INDEX[key]] += count
                depth += count
    pileup_tensor[BASE2INDEX['I1']] = max_ins_0
    pileup_tensor[BASE2INDEX['i1']] = max_ins_1
    pileup_tensor[BASE2INDEX['D1']] = max_del_0
    pileup_tensor[BASE2INDEX['d1']] = max_del_1
    denominator = depth if depth > 0 else 1
    pileup_list = sorted(list(pileup_dict.items()), key=lambda x: x[1], reverse=True)

    pass_snp_af = False
    pass_indel_af = False
    fast_mode = platform == 'ont' and fast_mode
    
    minimum_snp_af_for_candidate = minimum_snp_af_for_candidate if minimum_snp_af_for_candidate > 0 else param.threshold
    minimum_snp_af_for_candidate = max(minimum_snp_af_for_candidate, param.threshold_dict[platform]) if fast_mode else minimum_snp_af_for_candidate
    minimum_indel_af_for_candidate = minimum_indel_af_for_candidate if minimum_indel_af_for_candidate > 0 else param.threshold_dict[platform]

    # check whether first non reference candidate in the first position
    pass_af = len(pileup_list) and (pileup_list[0][0] != reference_base)
                                    
    for item, count in pileup_list:
        if item == reference_base:
            continue
        elif item[0] in 'ID':
            pass_indel_af = (pass_indel_af or (float(count) / denominator >= minimum_indel_af_for_candidate))
            continue
        if fast_mode:
            pass_snp_af = pass_snp_af or (float(count) / denominator >= minimum_snp_af_for_candidate and count >= 4)
        else:
            pass_snp_af = pass_snp_af or (float(count) / denominator >= minimum_snp_af_for_candidate)

    af = (float(pileup_list[1][1]) / denominator) if len(pileup_list) > 1 else 0.0
    af = (float(pileup_list[0][1]) / denominator) if len(pileup_list) >= 1 and pileup_list[0][
        0] != reference_base else af

    pileup_tensor[BASE2INDEX[reference_base]] = -1 * sum([pileup_tensor[BASE2INDEX[item]] for item in 'ACGT'])
    pileup_tensor[BASE2INDEX[reference_base.lower()]] = -1 * sum([pileup_tensor[BASE2INDEX[item]] for item in 'acgt'])

    pass_af = pass_af or pass_snp_af or pass_indel_af
    # add a return: base_counter for generating GVCF
    return pileup_tensor, alt_dict, af, depth, pass_af, pileup_list, max_del_length


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def CreateTensorPileup(args):
    """
    Create pileup tensor for pileup model training or calling.
    use slide window to scan the whole candidate regions, keep all candidates over specific minimum allelic frequency
    and minumum depth, use samtools mpileup to store pileup info for pileup tensor generation. Only scan candidate
    regions once, we could directly get all variant candidates.
    """
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    bed_file_path = args.bed_fn
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    samtools_execute_command = args.samtools
    bam_file_path = args.bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    tensor_can_output_path = args.tensor_can_fn
    is_bed_file_given = bed_file_path is not None
    minimum_af_for_candidate = args.threshold
    minimum_snp_af_for_candidate = args.snp_threshold
    minimum_indel_af_for_candidate = args.indel_threshold
    min_coverage = args.minCoverage
    platform = args.platform
    confident_bed_fn = args.confident_bed_fn
    is_confident_bed_file_given = confident_bed_fn is not None
    alt_fn = args.indel_fn
    extend_confident_bed_fn = args.extend_confident_bed_fn
    is_extend_confident_bed_file_given = extend_confident_bed_fn is not None
    min_mapping_quality = args.minMQ
    min_base_quality = args.minBQ
    fast_mode = args.fast_mode

    global test_pos
    test_pos = None
    # preparation for candidates near variants
    if not isfile("{}.fai".format(fasta_file_path)):
        sys.exit("Fasta index {}.fai doesn't exist.".format(fasta_file_path))

    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []
    tree, bed_start, bed_end = bed_tree_from(bed_file_path=extend_confident_bed_fn, contig_name=ctg_name,
                                             return_bed_region=True)

    fai_fn = file_path_from(fasta_file_path + ".fai", exit_on_not_found=True)
    if not is_bed_file_given and chunk_id is not None:
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

    if is_bed_file_given and chunk_id is not None:
        chunk_size = (bed_end - bed_start) // chunk_num + 1 if (bed_end - bed_start) % chunk_num else (
                                                                                                                  bed_end - bed_start) // chunk_num
        ctg_start = bed_start + 1 + chunk_size * chunk_id  # 0-base to 1-base
        ctg_end = ctg_start + chunk_size

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    if is_ctg_range_given:
        extend_start = ctg_start - no_of_positions
        extend_end = ctg_end + no_of_positions
        reads_regions.append(region_from(ctg_name=ctg_name, ctg_start=extend_start, ctg_end=extend_end))
        reference_start, reference_end = ctg_start - param.expandReferenceRegion, ctg_end + param.expandReferenceRegion
        reference_start = 1 if reference_start < 1 else reference_start
        ref_regions.append(region_from(ctg_name=ctg_name, ctg_start=reference_start, ctg_end=reference_end))
    elif is_ctg_name_given:
        reads_regions.append(region_from(ctg_name=ctg_name))
        ref_regions.append(region_from(ctg_name=ctg_name))
        reference_start = 1

    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools_execute_command,
        fasta_file_path=fasta_file_path,
        regions=ref_regions
    )

    if reference_sequence is None or len(reference_sequence) == 0:
        sys.exit("[ERROR] Failed to load reference seqeunce from file ({}).".format(fasta_file_path))

    if is_bed_file_given and ctg_name not in tree:
        sys.exit("[ERROR] ctg_name({}) not exists in bed file({}).".format(ctg_name, bed_file_path))

    # samtools mpileup options
    # reverse-del: deletetion in reverse strand marks as '#'
    min_mapping_quality = 0 if args.gvcf else min_mapping_quality
    min_base_quality = 0 if args.gvcf else min_base_quality
    max_depth = param.max_depth_dict[args.platform] if args.platform else args.max_depth
    mq_option = ' --min-MQ {}'.format(min_mapping_quality)
    bq_option = ' --min-BQ {}'.format(min_base_quality)
    flags_option = ' --excl-flags {}'.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    max_depth_option = ' --max-depth {}'.format(max_depth)
    bed_option = ' -l {}'.format(extend_confident_bed_fn) if is_extend_confident_bed_file_given else ""
    gvcf_option = ' -a' if args.gvcf else ""
    samtools_mpileup_process = subprocess_popen(
        shlex.split(
            "{} mpileup  {} -r {} --reverse-del".format(samtools_execute_command,
                                                        bam_file_path,
                                                        " ".join(reads_regions), )
            + mq_option + bq_option + bed_option + flags_option + max_depth_option + gvcf_option))

    if tensor_can_output_path != "PIPE":
        tensor_can_fpo = open(tensor_can_output_path, "wb")
        tensor_can_fp = subprocess_popen(shlex.split("{} -c".format(param.zstd)), stdin=PIPE, stdout=tensor_can_fpo)
    else:
        tensor_can_fp = TensorStdout(sys.stdout)

    # whether save all all alternative information, only for debug mode
    if alt_fn:
        alt_fp = open(alt_fn, 'w')

    pos_offset = 0
    pre_pos = -1
    tensor = [[]] * sliding_window_size
    candidate_position = []
    all_alt_dict = {}
    depth_dict = {}
    af_dict = {}

    # to generate gvcf, it is needed to record whole genome statistical information
    if args.gvcf:
        nonVariantCaller = variantInfoCalculator(gvcfWritePath=args.temp_file_dir, ref_path=args.ref_fn,
                                                 bp_resolution=args.bp_resolution, sample_name='.'.join(
                [args.sampleName, ctg_name, str(ctg_start), str(ctg_end)]), p_err=args.base_err,
                                                 gq_bin_size=args.gq_bin_size)

    confident_bed_tree = bed_tree_from(bed_file_path=confident_bed_fn, contig_name=ctg_name, bed_ctg_start=extend_start,
                                       bed_ctg_end=extend_end)

    for row in samtools_mpileup_process.stdout:
        columns = row.strip().split('\t')
        pos = int(columns[1])
        pileup_bases = columns[4]
        reference_base = reference_sequence[pos - reference_start].upper()
        valid_reference_flag = reference_base.upper() in ['A', 'T', 'C', 'G']
        if not valid_reference_flag and args.gvcf:
            nonVariantCaller.make_gvcf_online({}, push_current=True)
        if (ctg_start != None and ctg_end != None):
            within_flag = pos >= ctg_start and pos <= ctg_end
        elif (ctg_start != None and ctg_end == None):
            within_flag = pos >= ctg_start
        elif (ctg_start == None and ctg_end != None):
            within_flag = pos <= ctg_end
        else:
            within_flag = True
        if (args.gvcf and columns[3] == '0' and within_flag and valid_reference_flag):
            cur_site_info = {'chr': columns[0], 'pos': pos, 'ref': reference_base, 'n_total': 0, 'n_ref': 0}
            nonVariantCaller.make_gvcf_online(cur_site_info)
            continue

        # start with a new region, clear all sliding windows cache, avoid memory occupation
        if pre_pos + 1 != pos:
            pos_offset = 0
            tensor = [[]] * sliding_window_size
            candidate_position = []
        pre_pos = pos

        # a condition to skip some positions creating tensor,but return allele summary
        # allele count function
        pileup_tensor, alt_dict, af, depth, pass_af, pileup_list, max_del_length = generate_tensor(pos=pos,
                                                                                                   pileup_bases=pileup_bases,
                                                                                                   reference_sequence=reference_sequence,
                                                                                                   reference_start=reference_start,
                                                                                                   reference_base=reference_base,
                                                                                                   minimum_af_for_candidate=minimum_af_for_candidate,
                                                                                                   minimum_snp_af_for_candidate=minimum_snp_af_for_candidate,
                                                                                                   minimum_indel_af_for_candidate=minimum_indel_af_for_candidate,
                                                                                                   platform=platform,
                                                                                                   fast_mode=fast_mode)
        if args.gvcf and within_flag and valid_reference_flag:
            cur_n_total = 0
            cur_n_ref = 0
            for _key, _value in pileup_list:
                if (_key == reference_base):
                    cur_n_ref = _value
                cur_n_total += _value

            cur_site_info = {'chr': columns[0], 'pos': pos, 'ref': reference_base, 'n_total': cur_n_total,
                             'n_ref': cur_n_ref}
            nonVariantCaller.make_gvcf_online(cur_site_info)

        pass_confident_bed = not is_confident_bed_file_given or is_region_in(tree=confident_bed_tree,
                                                                             contig_name=ctg_name,
                                                                             region_start=pos - 1,
                                                                             region_end=pos + max_del_length + 1)  # 0-based
        if pass_confident_bed and reference_base in 'ACGT' and (pass_af and depth >= min_coverage):
            candidate_position.append(pos)
            all_alt_dict[pos] = alt_dict
            depth_dict[pos] = depth
            af_dict[pos] = af
        tensor[pos_offset] = pileup_tensor

        # save pileup tensor for each candidate position with nearby flanking_base_num bp distance
        pos_offset = (pos_offset + 1) % sliding_window_size
        if len(candidate_position) and pos - candidate_position[0] == flanking_base_num:
            center = candidate_position.pop(0)
            has_empty_tensor = sum([True for item in tensor if not len(item)])
            if not has_empty_tensor:
                depth = depth_dict[center]
                ref_seq = reference_sequence[center - (
                    flanking_base_num) - reference_start: center + flanking_base_num + 1 - reference_start]
                concat_tensor = tensor[pos_offset:] + tensor[0:pos_offset]

                alt_info = str(depth) + '-' + ' '.join(
                    [' '.join([item[0], str(item[1])]) for item in list(all_alt_dict[center].items())])
                l = "%s\t%d\t%s\t%s\t%s" % (
                    ctg_name,
                    center,
                    ref_seq,
                    " ".join(" ".join("%d" % x for x in innerlist) for innerlist in concat_tensor),
                    alt_info
                )
                tensor_can_fp.stdin.write(l)
                tensor_can_fp.stdin.write("\n")
                if alt_fn:
                    alt_info = ' '.join(
                        [' '.join([item[0], str(item[1])]) for item in list(all_alt_dict[center].items())])
                    alt_fp.write(
                        '\t'.join([ctg_name + ' ' + str(center), str(depth), alt_info, str(af_dict[center])]) + '\n')
                del all_alt_dict[center], depth_dict[center], af_dict[center]

    if (args.gvcf and len(nonVariantCaller.current_block) != 0):
        nonVariantCaller.write_to_gvcf_batch(nonVariantCaller.current_block, nonVariantCaller.cur_min_DP,
                                             nonVariantCaller.cur_raw_gq)

    samtools_mpileup_process.stdout.close()
    samtools_mpileup_process.wait()

    if tensor_can_output_path != "PIPE":
        tensor_can_fp.stdin.close()
        tensor_can_fp.wait()
        tensor_can_fpo.close()

    if alt_fn:
        alt_fp.close()


def main():
    parser = ArgumentParser(description="Generate 1-based variant candidates tensor using pileup alignments")

    parser.add_argument('--bam_fn', type=str, default="input.bam", required=True,
                        help="Sorted bam file input, default: %(default)s")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in provide candidate positions, bed regions work in intersection with each candidate positions with flanking 16 bp distance to speed up mpileup")

    parser.add_argument('--confident_bed_fn', type=str, default=None,
                        help="Call variant only in these regions, works in intersection with ctgName, ctgStart and ctgEnd, optional, default: as defined by ctgName, ctgStart and ctgEnd")

    parser.add_argument('--extend_confident_bed_fn', type=str, default=None,
                        help="Extended regions by confident bed regions to handle mpileup with candidates near provide bed regions, default extend 16 bp distance")

    parser.add_argument('--threshold', type=float, default=0.08,
                        help="Minimum allele frequence of the 1st non-reference allele for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--snp_threshold', type=float, default=0.0,
                        help="Minimum snp allele frequence of the 1st non-reference allele for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--indel_threshold', type=float, default=0.0,
                        help="Minimum indel allele frequence of the 1st non-reference allele for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--minCoverage', type=float, default=2,
                        help="Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="Minimum Base Quality. Base quality lower than the setting will be filtered, default: %(default)d")

    parser.add_argument('--ctgName', type=str, default="chr20",
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools verision >= 1.10 is required. default: %(default)s")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--indel_fn', type=str, default=None,
                        help="Output all alternative indel cigar, only use for analysis, default: %(default)s")

    parser.add_argument('--chunk_id', type=int, default=None,
                        help="Specific chunk id works with total chunk_num for parallel execution.")

    parser.add_argument('--chunk_num', type=int, default=None,
                        help="Total chunk number for parallel execution. Each chunk refer to a smaller reference regions.")

    parser.add_argument('--max_depth', type=int, default=144,
                        help="Maximum pileup depth to be processed. default: %(default)s")

    parser.add_argument('--zstd', type=str, default=param.zstd,
                        help="Path to the 'zstd' compression, default: %(default)s")

    parser.add_argument('--test_pos', type=int, default=0,
                        help="Test in specific candidate position. Only use for analysis, deprecated")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select specific platform for variant calling. Optional: 'ont,pb,illumina', default: %(default)s")

    parser.add_argument('--fast_mode', type=str2bool, default=False,
                        help="Ignore low allelic frequency <= 0.15 snp calling for ont platform,  optional")

    # gvcf options
    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Whether to generate gvcf")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='Default estimated base error rate')

    parser.add_argument('--gq_bin_size', default=5, type=int,
                        help='Default gq bin size for merge non-variant block')

    parser.add_argument('--temp_file_dir', type=str, default="./",
                        help="Temporary directory for save the temporary gvcf files")

    parser.add_argument('--bp_resolution', action='store_true',
                        help="Whether bp resolution for GVCF")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    CreateTensorPileup(args)


if __name__ == "__main__":
    main()
