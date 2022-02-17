import sys
import shlex
import os
import re
import json
import logging
import random
from subprocess import PIPE
from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict, OrderedDict

import shared.param_f as param
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool, vcf_candidates_from
from shared.interval_tree import bed_tree_from, is_region_in


logging.basicConfig(format='%(message)s', level=logging.INFO)
BASES = set(list(BASE2NUM.keys()) + ["-"])
no_of_positions = param.no_of_positions
flanking_base_num = param.flankingBaseNum
channel_size = param.channel_size
BASE2NUMBER = dict(zip("ACGTURYSWKMBDHVN-", (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 4)))
NORMALIZE_NUM = param.NORMALIZE_NUM
MAX_BQ = 40.0
MAX_MQ = 60.0
MAX_AF = 1.0
STRAND_0 = 100
STRAND_1 = 50
HAP_TYPE = dict(zip((1, 0, 2), (30, 60, 90)))  # hap1 UNKNOWN H2
ACGT_NUM = dict(zip("ACGT+-*#N", (100, 25, 75, 50, -50, -100, 0, 0, 100)))
cigarRe = r"(\d+)([MIDNSHP=X])"


def get_cigar_tuple(cigar):
    cigar_tuple = []
    for m in re.finditer(cigarRe, cigar):
        cigar_tuple.append((m.group(2), (int(m.group(1)))))

    return cigar_tuple


def _normalize_bq(x):
    return int(NORMALIZE_NUM * min(x, MAX_BQ) / MAX_BQ)


def _normalize_mq(x):
    return int(NORMALIZE_NUM * min(x, MAX_MQ) / MAX_MQ)


def _normalize_af(x):
    return int(NORMALIZE_NUM * min(x, MAX_AF) / MAX_AF)


def edit_distance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

try:
    #pypy3 -m pip install python-Levenshtein==0.12.2
    from Levenshtein import distance as edit_distance
    distance = edit_distance
except:
    # python version is much slower
    distance = edit_distance


class Alignment(object):
    def __init__(self, cigartuples, query_sequence, reference_start, read_name):
        self.cigartuples = cigartuples
        self.query_sequence = query_sequence
        self.reference_start = reference_start
        self.read_name = read_name


class Variant(object):
    def __init__(self, position, reference_allele=None, alternative_allele=None, genotype=None, phase_set=None):
        self.position = position
        self.reference_allele = reference_allele
        self.alternative_allele = alternative_allele
        self.genotype = genotype
        self.phase_set = phase_set


class Read(object):
    def __init__(self, read_start, strand, MQ, ):
        # self.pos_channel = defaultdict(list)
        self.read_start = read_start
        self.pos_channel = defaultdict(list)
        self.pos_alt = defaultdict(str)
        self.pos_ins = defaultdict(str)
        self.read_end = None
        self.strand = STRAND_1 if strand == True else STRAND_0
        self.MQ = _normalize_mq(MQ)
        self.HP = 60


class Position(object):
    def __init__(self):

        self.read_name_list = []
        self.base_list = []
        self.raw_base_quality = []
        self.raw_mapping_quality = []
        self.depth = 0
        self.read_channel = None
        self.mapping_quality = None
        self.update_info = False
        self.read_info = defaultdict()
        self.alt_dict = defaultdict(int)


def _iterate_cigar(variants, j, bam_read, cigartuples):
    """
    Iterate over the CIGAR of the given bam_read and variants[j:] in lockstep.
    Yield tuples (index, i, consumed, query_pos) where index is into the variants list
    i and consumed describe the split position in the cigar
    bam_read -- a pysam.AlignedSegment
    variants -- list of variants (VcfVariant objects)
    j -- index of the first variant (in the variants list) to check
    """
	# Skip variants that are located to the left of the read
    ref_pos = bam_read.reference_start
    n = len(variants)
    query_pos = 0
    while j < n and variants[j].position < ref_pos:
        j += 1

    # Iterate over the CIGAR sequence (defining the alignment) and variant list in lockstep
    for i, (cigar_op, length) in enumerate(cigartuples):
        if j < n:
            v_position = variants[j].position
        if cigar_op in "MX=":  # M, X, = operators (match)
            # Iterate over all variants that are in this matching region
            while j < n and v_position < ref_pos + length:
                assert v_position >= ref_pos
                yield (j, i, v_position - ref_pos, query_pos + v_position - ref_pos)
                j += 1
                if j < n:
                    v_position = variants[j].position
            query_pos += length
            ref_pos += length
        elif cigar_op == 'I':  # I operator (insertion)
            # TODO it should work to *not* handle the variant here, but at the next M or D region
            if j < n and v_position == ref_pos:
                yield (j, i, 0, query_pos)
                j += 1
                if j < n:
                    v_position = variants[j].position
            query_pos += length
        elif cigar_op == 'D':  # D operator (deletion)
            # Iterate over all variants that are in this deleted region
            while j < n and v_position < ref_pos + length:
                assert v_position >= ref_pos
                yield (j, i, v_position - ref_pos, query_pos)
                j += 1
                if j < n:
                    v_position = variants[j].position
            ref_pos += length
        elif cigar_op == 'N':  # N operator (reference skip)
            # Iterate over all variants that are in this skipped region
            while j < n and v_position < ref_pos + length:
                assert v_position >= ref_pos
                j += 1
                if j < n:
                    v_position = variants[j].position
            ref_pos += length
        elif cigar_op == 'S':  # S operator (soft clipping)
            query_pos += length
        elif cigar_op in 'HP':  # H or P (hard clipping or padding)
            pass
        else:
            raise ValueError("Unsupported CIGAR operation: {}".format(cigar_op))


def split_cigar(cigar, i, consumed):
    """
    Split a CIGAR into two parts. i and consumed describe the split position.
    i is the element of the cigar list that should be split, and consumed says
    at how many operations to split within that element.

    The CIGAR is given as a list of (operation, length) pairs.

    i -- split at this index in cigar list
    consumed -- how many cigar ops at cigar[i] are to the *left* of the
        split position

    Return a tuple (left, right).

    Example:
    Assume the cigar is 3M 1D 6M 2I 4M.
    With i == 2 and consumed == 5, the cigar is split into
    3M 1D 5M and 1M 2I 4M.
    """
    middle_op, middle_length = cigar[i]
    assert consumed <= middle_length
    if consumed > 0:
        left = cigar[:i] + [(middle_op, consumed)]
    else:
        left = cigar[:i]
    if consumed < middle_length:
        right = [(middle_op, middle_length - consumed)] + cigar[i + 1 :]
    else:
        right = cigar[i + 1 :]
    return left, right

def cigar_prefix_length(cigar, reference_bases):
    """
    Given a prefix of length reference_bases relative to the reference, how
    long is the prefix of the read? In other words: If reference_bases on
    the reference are consumed, how many bases on the query does that
    correspond to?

    If the position is within or at the end of an insertion (which do not
    consume bases on the reference), then the number of bases up to the
    beginning of the insertion is reported.

    Return a pair (reference_bases, query_bases) where the value for
    reference_bases may be smaller than the requested one if the CIGAR does
    not cover enough reference bases.

    Reference skips (N operators) are treated as the end of the read. That
    is, no positions beyond a reference skip are reported.
    """
    ref_pos = 0
    query_pos = 0
    for op, length in cigar:
        if op in "MX=":  # M, X, =
            ref_pos += length
            query_pos += length
            if ref_pos >= reference_bases:
                return (reference_bases, query_pos + reference_bases - ref_pos)
        elif op == 'D':  # D
            ref_pos += length
            if ref_pos >= reference_bases:
                return (reference_bases, query_pos)
        elif op == 'I':  # I
            query_pos += length
        elif op in "SH":  # soft or hard clipping
            pass
        elif op == 'N':  # N
            # Always stop at reference skips
            return (reference_bases, query_pos)
        else:
            assert False, "unknown CIGAR operator"
    assert ref_pos < reference_bases
    return (ref_pos, query_pos)


def realign(
    variant,
    bam_read,
    cigartuples,
    i,
    consumed,
    query_pos,
    reference,
    overhang,
    reference_start
):
    """
    Realign a read to the two alleles of a single variant.
    i and consumed describe where to split the cigar into a part before the
    variant position and into a part starting at the variant position, see split_cigar().

    variant -- VcfVariant
    bam_read -- the AlignedSegment
    cigartuples -- the AlignedSegment.cigartuples property (accessing it is expensive, so re-use it)
    i, consumed -- see split_cigar method
    query_pos -- index of the query base that is at the variant position
    reference -- the reference as a str-like object (full chromosome)
    overhang -- extend alignment by this many bases to left and right
    use_affine -- if true, use affine gap costs for realignment
    gap_start, gap_extend -- if affine_gap=true, use these parameters for affine gap cost alignment
    default_mismatch -- if affine_gap=true, use this as mismatch cost in case no base qualities are in bam
    """
    # Do not process symbolic alleles like <DEL>, <DUP>, etc.
    if variant.alternative_allele.startswith("<"):
        return None, None

    left_cigar, right_cigar = split_cigar(cigartuples, i, consumed)

    left_ref_bases, left_query_bases = cigar_prefix_length(
        left_cigar[::-1], overhang
    )
    right_ref_bases, right_query_bases = cigar_prefix_length(
        right_cigar, len(variant.reference_allele) + overhang
    )

    assert variant.position - left_ref_bases >= 0
    # assert variant.position + right_ref_bases <= len(reference)

    query = bam_read.query_sequence[
        query_pos - left_query_bases : query_pos + right_query_bases
    ]
    ref = reference[variant.position - left_ref_bases - reference_start + 1: variant.position + right_ref_bases - reference_start + 1]
    alt = (
        reference[variant.position - left_ref_bases - reference_start + 1: variant.position - reference_start + 1]
        + variant.alternative_allele
        + reference[
            variant.position
            + len(variant.reference_allele) - reference_start + 1: variant.position
            + right_ref_bases - reference_start + 1
        ]
    )

    base_qual_score = 30
    distance_ref = distance(query, ref)
    distance_alt = distance(query, alt)

    if distance_ref < distance_alt:
        return 0, base_qual_score  # detected REF
    elif distance_ref > distance_alt:
        return 1, base_qual_score  # detected ALT
    else:
        return None, None  # cannot decide


def detect_alleles_by_alignment(
    variants,
    j,
    bam_read,
    reference,
    reference_start=1,
    overhang=10,
):
    """
    Detect which alleles the given bam_read covers. Detect the correct
    alleles of the variants that are covered by the given bam_read.

    Yield tuples (position, allele, quality).

    variants -- list of variants (VcfVariant objects)
    j -- index of the first variant (in the variants list) to check
    """
    # Accessing bam_read.cigartuples is expensive, do it only once
    cigartuples = bam_read.cigartuples

    # For the same reason, the following check is here instad of
    # in the _usable_alignments method
    if not cigartuples:
        return

    for index, i, consumed, query_pos in _iterate_cigar(variants, j, bam_read, cigartuples):
        allele, quality = realign(
            variants[index],
            bam_read,
            cigartuples,
            i,
            consumed,
            query_pos,
            reference,
            overhang,
            reference_start
        )
        if allele in (0, 1):
            yield (index, allele, quality)

def phredscore2raw_score(qual):
    return ord(qual) - 33

def raw_score2phredscore(s):
    return chr(s+33)

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


def update_read_info(read_info, pos, base_info, ref_base):
    base, indel, bq = base_info
    ins_base = ""
    if pos not in read_info.pos_channel:
        read_info.pos_channel[pos] = [0] * channel_size
    if base == "#" or ref_base not in 'ACGTN':
        return
    REF_BASE = ACGT_NUM[ref_base]
    ALT_BASE = 0
    if indel != '':
        ALT_BASE = ACGT_NUM[indel[0]]
        if indel[0] == '+':
            ins_base = indel[1:]
    elif (base != ref_base and base in 'ACGT'):
        base = evc_base_from(base)
        ALT_BASE = ACGT_NUM[base]
    bq = 0 if bq == '' else _normalize_bq(phredscore2raw_score(bq))
    read_info.pos_channel[pos][:5] = REF_BASE, ALT_BASE, read_info.strand, read_info.MQ, bq
    read_info.pos_channel[pos][7] = read_info.HP
    if ins_base != "":
        read_info.pos_ins[pos] = ins_base


def sorted_by_hap_read_name(center_pos, all_read_dict, hap_dict, platform):
    """
    Sort by reads haplotype after haplotag reads otherwise sort by read start position.
    center_pos: define the center candidate position for proccessing.
    haplotag_dict: dictionary (read name : hap type) which keep the read name and haplotype mapping.
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position .
    hap_dict: similar to haplotag_dict, dictionary (pos: pos info) which keep the read name and haplotype mapping,
    while haplotype information directly acquire from BAM HP tag.
    platform: select maximum depth for each platform.
    """
    all_nearby_read_name = []
    start_pos, end_pos = center_pos - flanking_base_num, center_pos + flanking_base_num
    for read_name, read in all_read_dict.items():
        read_start, read_end = read.read_start, read.read_end
        if read_start > end_pos:
            break
        if read_end <= start_pos:
            continue
        all_nearby_read_name.append(read_name)

    all_nearby_read_name = list(OrderedDict.fromkeys(all_nearby_read_name))  # have sorted by order
    matrix_depth = param.matrix_depth_dict[platform]
    if len(all_nearby_read_name) > matrix_depth:
        # set same seed for reproducibility
        random.seed(0)
        indices = random.sample(range(len(all_nearby_read_name)), matrix_depth)
        all_nearby_read_name = [all_nearby_read_name[i] for i in sorted(indices)]
    sorted_read_name_list = []
    for order, read_name in enumerate(all_nearby_read_name):
        hap = hap_dict[read_name]  # no phasing is 0
        sorted_read_name_list.append((hap, order, read_name))

    sorted_read_name_list = sorted(sorted_read_name_list)
    return sorted_read_name_list

def get_alt_info(center_pos, pileup_dict, ref_seq, reference_sequence, reference_start, hap_dict):
    """
    Get alternative information for representation unification, keep all read level alignment information including phasing info.
    center_pos: center position for processing, default window size = no_of_positions = flankingBaseNum + 1 + flankingBaseNum
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position .
    ref_seq: chunked reference sequence in window, start: center pos - flankingBaseNum, end: center + flankingBaseNum + 1.
    reference_sequence: reference sequence index by contig:start-end. 0-based.
    reference_base: upper reference base for cigar calculation.
    reference_start: upper reference base for cigar calculation.
    hap_dict: dictionary (pos: pos info) which keep the read name and haplotype mapping.
    """

    reference_base = ref_seq[flanking_base_num]
    alt_read_name_dict = defaultdict(set)
    depth = 0
    for (base, indel), read_name in zip(pileup_dict[center_pos].base_list, pileup_dict[center_pos].read_name_list):
        if base in "#*":
            alt_read_name_dict['*'].add(read_name)
            depth += 1
            continue
        depth += 1
        if base.upper() == reference_base and indel == '':
            alt_read_name_dict['R'].add(read_name)
        if indel != '':
            if indel[0] == '+':
                indel = 'I' + base.upper() + indel.upper()[1:]
            else:
                del_bases_num = len(indel[1:])
                del_ref_bases = reference_sequence[
                                center_pos - reference_start + 1:center_pos - reference_start + del_bases_num + 1]
                indel = 'D' + del_ref_bases
            alt_read_name_dict[indel].add(read_name)

        if indel == '' and base.upper() != reference_base:
            alt_read_name_dict['X' + base.upper()].add(read_name)

    for alt_type, read_name_set in list(alt_read_name_dict.items()):
        alt_read_name_dict[alt_type] = ' '.join(
            [read_name + '_' + str(hap_dict[read_name]) for read_name in list(read_name_set)])

    alt_info = str(depth) + '\t' + json.dumps(alt_read_name_dict)

    return alt_info


def generate_tensor(ctg_name, center_pos, sorted_read_name_list, pileup_dict, ref_seq, reference_sequence,
                    reference_start, platform, confident_bed_tree, add_no_phasing_data_training, no_phasing_data_training_proportion=1.0, all_read_dict=None):
    """
    Generate full alignment input tensor
    ctg_name: provided contig name.
    center_pos: center position for full alignment generation, default window size = no_of_positions =
    flankingBaseNum + 1 + flankingBaseNum
    sorted_read_name_list: read name list which have been sorted by read start position and haplotype.
    pileup_dict: dictionary (pos: pos info) which keep read information that cover specific position .
    ref_seq: chunked reference sequence in window, start: center pos - flankingBaseNum, end: center + flankingBaseNum + 1.
    reference_sequence: reference sequence index by contig:start-end. 0-based.
    reference_base: upper reference base for cigar calculation.
    reference_start: upper reference base for cigar calculation.
    platform: platform for tensor shape, ont give a larger maximum depth compared with pb and illumina.
    confident_bed_tree: dictionary (contig name : intervaltree) for fast region query.
    add_no_phasing_data_training: boolean option to decide whether add no phasing data in training, we will
    resort the read and remove haplotype info when using this option.
    """

    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    reference_base = ref_seq[flanking_base_num]
    tensor_depth = len(sorted_read_name_list)
    if tensor_depth == 0:
        return None, None
    tensor = [[[0] * tensor_shape[2] for _ in range(tensor_shape[1])] for _ in range(tensor_depth)]
    start_pos, end_pos = center_pos - flanking_base_num, center_pos + flanking_base_num + 1
    insert_tuple = []
    depth = max(1, pileup_dict[center_pos].depth)
    alt_dict = pileup_dict[center_pos].alt_dict

    pass_confident_bed = not len(confident_bed_tree) or is_region_in(confident_bed_tree, ctg_name,
                                                                     center_pos - 2,
                                                                     center_pos + 1)
    if not pass_confident_bed:
        return None, None

    for p in range(start_pos, end_pos):

        for read_idx, read_name_info in enumerate(sorted_read_name_list):
            hap, _, read_name = read_name_info
            offset = p - start_pos
            pos_channel_dict = all_read_dict[read_name].pos_channel
            if p in pos_channel_dict:
                tensor[read_idx][offset] = pos_channel_dict[p]
                ins_base = all_read_dict[read_name].pos_ins[p]
                if ins_base != '' and p < end_pos - 1:
                    insert_tuple.append((read_idx, offset, ins_base))

    for read_idx, offset, ins_base in insert_tuple:
        for ins_idx in range(min(len(ins_base), no_of_positions - offset)):
            tensor[read_idx][ins_idx + offset][6] = ACGT_NUM[ins_base[ins_idx]]

    for row_idx, (hap, _, read_name) in enumerate(sorted_read_name_list):
        af_num = 0
        if center_pos in all_read_dict[read_name].pos_alt:
            base, indel = all_read_dict[read_name].pos_alt[center_pos]
            base_upper = base.upper()
            if indel != '':
                if indel[0] == '+':
                    af_num = alt_dict[indel] / depth if indel in alt_dict else af_num
                else:
                    af_num = alt_dict[indel] / depth if indel in alt_dict else af_num
            elif base in alt_dict:
                af_num = alt_dict[base_upper] / depth
        af_num = _normalize_af(af_num) if af_num != 0 else af_num

        for p in range(no_of_positions):
            if tensor[row_idx][p][2] != 0:  # skip all del #*
                tensor[row_idx][p][5] = af_num


    alt_info = []
    for alt_type, alt_count in alt_dict.items():
        if alt_type[0] == '+':
            alt_info.append(['I' + reference_base + alt_type[1:].upper(), str(alt_count)])
        elif alt_type[0] == '-':
            del_bases_num = len(alt_type[1:])
            del_ref_bases = reference_sequence[
                            center_pos - reference_start + 1:center_pos - reference_start + del_bases_num + 1]
            alt_info.append(['D' + del_ref_bases, str(alt_count)])
        else:
            alt_info.append(['X' + alt_type, str(alt_count)])

    alt_info = str(depth) + '-' + ' '.join([' '.join([item[0], str(item[1])]) for item in alt_info])
    tensor_string_list = [" ".join((" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in tensor)]

    if add_no_phasing_data_training:
        all_hap = [item[0] for item in sorted_read_name_list]
        # skip if no phased reads exist
        if sum(all_hap) != 0:
            random.seed(int(center_pos))
            if random.random() <= no_phasing_data_training_proportion:
                raw_read_name_index_mapping = [item[1] for item in sorted(
                    [(item[1], read_idx) for read_idx, item in enumerate(sorted_read_name_list)])]
                no_phasing_tensor = [tensor[read_idx] for read_idx in raw_read_name_index_mapping]
                for row_idx in range(len(no_phasing_tensor)):
                    for p in range(no_of_positions):
                        if tensor[row_idx][p][7] > 0:
                            tensor[row_idx][p][7] = HAP_TYPE[0]

                no_phasing_tensor_string = " ".join(
                    (" ".join(" ".join(str(x) for x in innerlist) for innerlist in outerlist)) for outerlist in
                    no_phasing_tensor)
                tensor_string_list.append(no_phasing_tensor_string)
    return '\n'.join(["%s\t%d\t%s\t%s\t%s" % (
        ctg_name,
        center_pos,
        ref_seq,
        tensor_string,
        alt_info
    ) for tensor_string in tensor_string_list]), alt_info


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def CreateTensorFullAlignment(args):
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    full_aln_regions = args.full_aln_regions
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    need_haplotaging = args.need_haplotaging
    samtools_execute_command = args.samtools
    bam_file_path = args.bam_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    tensor_can_output_path = args.tensor_can_fn
    is_full_aln_regions_given = full_aln_regions is not None
    phasing_info_in_bam = args.phasing_info_in_bam
    phasing_window_size = args.phasing_window_size
    extend_bp = param.extend_bp
    unify_repre = args.unify_repre
    minimum_af_for_candidate = args.min_af
    minimum_snp_af_for_candidate = args.snp_min_af
    minimum_indel_af_for_candidate = args.indel_min_af
    min_coverage = args.minCoverage
    platform = args.platform
    confident_bed_fn = args.bed_fn
    is_confident_bed_file_given = confident_bed_fn is not None
    phased_vcf_fn = args.phased_vcf_fn
    alt_fn = args.indel_fn
    extend_bed = args.extend_bed
    is_extend_bed_file_given = extend_bed is not None
    min_mapping_quality = args.minMQ
    min_base_quality = args.minBQ
    unify_repre_fn = args.unify_repre_fn
    add_no_phasing_data_training = args.add_no_phasing_data_training
    vcf_fn = args.vcf_fn
    is_known_vcf_file_provided = vcf_fn is not None
    no_phasing_data_training_proportion = args.no_phasing_data_training_proportion
    global test_pos
    test_pos = None

    candidates_set = set()
    variants = []
    pileup_dict = defaultdict(str)
    need_process_pos_set = set()
    add_read_regions = True
    if full_aln_regions:

        """
        If given full alignment bed regions, all candidate positions will be directly selected from each row, define as 
        'ctg start end', where 0-based center position is the candidate for full alignment calling.
        if 'need_haplotaging' option enables, full alignment bed regions will also include nearby heterozygous snp candidates for reads
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

            if platform == "ilmn":
                continue
            if len(row) > 3:  # hete snp positions
                center_pos = position + extend_bp + 1
                ref_base, alt_base, genotype, phase_set = row[3].split('-')
                variant = Variant(position=center_pos-1, reference_allele=ref_base, alternative_allele=alt_base,
                                                         genotype=int(genotype), phase_set=phase_set)
                if phased_vcf_fn is None:
                    variants.append(variant)
            else:
                center = position + (end - position) // 2 - 1
                candidates_set.add(center)
                for p in list(range(center - flanking_base_num, center + flanking_base_num + 1)):
                    need_process_pos_set.add(p)
                    pileup_dict[p] = Position()

        candidate_file_path_output.close()
        candidate_file_path_process.wait()

    if platform == 'ilmn' and bam_file_path == "PIPE":
        add_read_regions = False

    fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')

    if is_known_vcf_file_provided:
        known_variants_list = vcf_candidates_from(vcf_fn=vcf_fn, contig_name=ctg_name)
        candidates_set = set(known_variants_list)
    if not full_aln_regions and chunk_id is not None:

        """
        Whole genome calling option, acquire contig start end position from reference fasta index(.fai), then split the
        reference accroding to chunk id and total chunk numbers.
        """
        if is_confident_bed_file_given:
            # consistent with pileup generation, faster to extract tensor using bed region
            tree, bed_start, bed_end = bed_tree_from(bed_file_path=extend_bed,
                                                     contig_name=ctg_name,
                                                     return_bed_region=True)

            chunk_size = (bed_end - bed_start) // chunk_num + 1 if (bed_end - bed_start) % chunk_num else (
                                                                                                                      bed_end - bed_start) // chunk_num
            ctg_start = bed_start + 1 + chunk_size * chunk_id  # 0-base to 1-base
            ctg_end = ctg_start + chunk_size
        else:
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

        # for illumina platform, the reads alignment is acquired after reads realignment from ReadsRealign.py
        if platform == 'ilmn' and bam_file_path != "PIPE":
            bam_file_path += '.{}_{}'.format(ctg_start, ctg_end)
            add_read_regions = False
        if bam_file_path == "PIPE":
            add_read_regions = False

    if need_haplotaging and phased_vcf_fn and os.path.exists(phased_vcf_fn):
        # if need_phasing option enables, scan the phased vcf file and store the heterozygous SNP candidates from each phase set
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

            variant = Variant(position=pos - 1, reference_allele=ref_base, alternative_allele=alt_base,
                           genotype=int(genotype), phase_set=phase_set)
            variants.append(variant)


    # preparation for candidates near variants
    candidates_set = set([item for item in candidates_set if item >= ctg_start and item <= ctg_end])
    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []

    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    extend_start, extend_end = None, None
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
        sys.exit("[ERROR] Failed to load reference sequence from file ({}).".format(fasta_file_path))

    stdin = None if bam_file_path != "PIPE" else sys.stdin
    bam_file_path = bam_file_path if bam_file_path != "PIPE" else "-"

    mq_option = ' -q {}'.format(min_mapping_quality) if min_mapping_quality > 0 else ""
    flags_option = ' -F {}'.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    bed_option = ' -L {}'.format(
        extend_bed) if is_extend_bed_file_given and platform != 'ilmn' else ""
    bed_option = ' -L {}'.format(full_aln_regions) if is_full_aln_regions_given and platform != 'ilmn' else bed_option
    reads_regions_option = ' ' + " ".join(reads_regions) if add_read_regions else ""

    samtools_command = "{} view {} ".format(samtools_execute_command, bam_file_path) + \
                         reads_regions_option + mq_option + bed_option + flags_option

    samtools_view_process = subprocess_popen(
        shlex.split(samtools_command), stdin=stdin)

    if not unify_repre:
        if tensor_can_output_path != "PIPE":
            tensor_can_fpo = open(tensor_can_output_path, "wb")
            tensor_can_fp = subprocess_popen(shlex.split("{} -c".format(args.zstd)), stdin=PIPE, stdout=tensor_can_fpo)
        else:
            tensor_can_fp = TensorStdout(sys.stdout)
    else:
        if unify_repre_fn != "PIPE":
            label_fp = open(unify_repre_fn, 'w')
        else:
            label_fp = sys.stdout
    if alt_fn:
        output_alt_fn = alt_fn
        alt_fp = open(output_alt_fn, 'w')

    hap_dict = defaultdict(int)
    all_read_dict = OrderedDict()
    extend_bp_distance = no_of_positions + param.extend_bp
    confident_bed_tree = bed_tree_from(bed_file_path=confident_bed_fn,
                                       contig_name=ctg_name,
                                       bed_ctg_start=extend_start,
                                       bed_ctg_end=extend_end)

    for p in pileup_dict:
        pileup_dict[p].ref_base = reference_sequence[p - reference_start]

    def samtools_view_generator_from(samtools_view_process):
        candidates_list = sorted(list(candidates_set))
        current_pos_index = 0
        varaint_current_pos = 0  # index of variants
        for row_id, row in enumerate(samtools_view_process.stdout):

            if row[0] == '@':
                continue
            columns = row.strip().split()
            read_name = columns[0]
            FLAG = int(columns[1])
            POS = int(columns[3])
            MAPQ = int(columns[4])
            CIGAR = columns[5]
            SEQ = columns[9].upper()
            BQ = columns[10]
            reference_position = POS - 1
            STRAND = (16 == (FLAG & 16))
            query_pos = 0

            read_info = Read(read_start=POS, strand=STRAND, MQ=MAPQ)
            cigartuples = get_cigar_tuple(CIGAR)

            if need_haplotaging:
                # haplotaging alignment here
                # 20 is minMQ in whatshap haplotaging
                if MAPQ < 20:
                    haplotype = 0
                else:

                    bam_read = Alignment(reference_start=reference_position,cigartuples=cigartuples, query_sequence=SEQ, read_name=read_name)
                    while (
                            varaint_current_pos < len(variants)
                            and variants[varaint_current_pos].position < reference_position
                    ):
                        varaint_current_pos += 1
                    detected = detect_alleles_by_alignment(
                        variants,
                        varaint_current_pos,
                        bam_read,
                        reference_sequence,
                        reference_start
                    )
                    haplotype_costs = defaultdict(int)
                    for j, allele, quality in detected:
                        phase_set = variants[j].phase_set
                        hete_hap_type = variants[j].genotype
                        if allele + 1 == hete_hap_type:
                            haplotype_costs[phase_set] += 1
                        else:
                            haplotype_costs[phase_set] -= 1

                    best_haplotype_cost = sorted(list(haplotype_costs.items()), key=lambda x: -abs(x[1]))
                    if len(best_haplotype_cost) == 0 or best_haplotype_cost[0][1] == 0:
                        haplotype = 0
                    else:
                        phaseset, quality = best_haplotype_cost[0]
                        haplotype = 1 if quality > 0 else 2

                hap_dict[read_name] = haplotype
                if haplotype > 0:
                    read_info.HP = HAP_TYPE[haplotype]

            base_list = []
            pos_list = []
            ref_pos = POS
            for i, (cigar_op, length) in enumerate(cigartuples):
                if cigar_op in 'MX=':
                    for p in range(ref_pos, ref_pos + length):
                        if p in need_process_pos_set:
                            base_list.append([SEQ[query_pos], "", BQ[query_pos]])
                            pos_list.append(p)
                        query_pos += 1
                    ref_pos += length
                elif cigar_op == 'D':
                    if p in need_process_pos_set:
                        base_list[-1][1] = '-' + 'N' * length
                    for p in range(ref_pos, ref_pos + length):
                        if p in need_process_pos_set:
                            base_list.append(["#", "", ""])
                            pos_list.append(p)
                    ref_pos += length
                elif cigar_op == 'I':
                    if p in need_process_pos_set:
                        base_list[-1][1] = '+' + SEQ[query_pos:query_pos+length]
                    query_pos += length
                elif cigar_op == 'N':
                    ref_pos += length
                elif cigar_op == 'S':
                    query_pos += length
                elif cigar_op in "HP":
                    pass
                else:
                    raise ValueError("Unsupported CIGAR operation: {}".format(cigar_op))

            read_info.read_end = ref_pos

            for p, base_info in zip(pos_list, base_list):
                ref_base = reference_sequence[p - reference_start]
                if p in candidates_set:
                    if base_info[0] != ref_base and (base_info[0] != "#") or base_info[1] != "":
                        if base_info[1] == "":
                            pileup_dict[p].alt_dict[base_info[0]] += 1
                        else:
                            pileup_dict[p].alt_dict[base_info[1]] += 1
                    read_info.pos_alt[p] = base_info[:2]

                update_read_info(read_info=read_info, pos=p, base_info=base_info, ref_base=ref_base.upper())
                pileup_dict[p].depth += 1
            all_read_dict[read_name] = read_info

            if current_pos_index < len(candidates_list) and POS - candidates_list[
                current_pos_index] > extend_bp_distance:
                yield candidates_list[current_pos_index]

                current_pos_index += 1
        while current_pos_index != len(candidates_list):
            yield candidates_list[current_pos_index]
            current_pos_index += 1

    samtools_pileup_generator = samtools_view_generator_from(samtools_view_process)

    for p_id, pos in enumerate(samtools_pileup_generator):
        if pos not in pileup_dict:
            continue
        if p_id > 0 and p_id % 100 == 0:
            for read_name in list(all_read_dict.keys()):
                if all_read_dict[read_name].read_end < pos - extend_bp_distance:
                    del all_read_dict[read_name]
            for pre_pos in sorted(pileup_dict.keys()):
                if pos - pre_pos > extend_bp_distance:
                    del pileup_dict[pre_pos]
                else:
                    break

        sorted_read_name_list = sorted_by_hap_read_name(pos, all_read_dict, hap_dict, platform)
        ref_seq = reference_sequence[
                  pos - reference_start - flanking_base_num: pos - reference_start + flanking_base_num + 1].upper()

        if not unify_repre:
            tensor, alt_info = generate_tensor(ctg_name=ctg_name,
                                               center_pos=pos,
                                               sorted_read_name_list=sorted_read_name_list,
                                               pileup_dict=pileup_dict,
                                               ref_seq=ref_seq,
                                               reference_sequence=reference_sequence,
                                               reference_start=reference_start,
                                               platform=platform,
                                               confident_bed_tree=confident_bed_tree,
                                               add_no_phasing_data_training=add_no_phasing_data_training,
                                               no_phasing_data_training_proportion=no_phasing_data_training_proportion,
                                               all_read_dict=all_read_dict)


            if not tensor:
                continue

            tensor_can_fp.stdin.write(tensor)
            tensor_can_fp.stdin.write("\n")
            if alt_fn:
                alt_info = alt_info.replace('-', '\t')
                alt_fp.write('\t'.join([ctg_name + ' ' + str(pos), alt_info]) + '\n')

        if unify_repre:
            label_info = get_alt_info(center_pos=pos,
                                      pileup_dict=pileup_dict,
                                      ref_seq=ref_seq,
                                      reference_sequence=reference_sequence,
                                      reference_start=reference_start,
                                      hap_dict=hap_dict)
            label_fp.write('\t'.join([ctg_name + ' ' + str(pos), label_info]) + '\n')

    samtools_view_process.stdout.close()
    samtools_view_process.wait()

    if not unify_repre and tensor_can_output_path != "PIPE":
        tensor_can_fp.stdin.close()
        tensor_can_fp.wait()
        tensor_can_fpo.close()

    if alt_fn:
        alt_fp.close()

    if unify_repre and unify_repre_fn != "PIPE":
        label_fp.close()


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
    parser.add_argument('--minCoverage', type=float, default=param.min_coverage,
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
    parser.add_argument('--need_haplotaging', action='store_true',
                        help=SUPPRESS)

    ## Apply read realignment for illumina platform. Greatly boost indel performance in trade of running time
    parser.add_argument('--need_realignment', action='store_true',
                        help=SUPPRESS)

    ## No phasing data training proportion in phased input
    parser.add_argument('--no_phasing_data_training_proportion', default=param.no_phasing_data_training_proportion, type=float,
                        help=SUPPRESS)

    parser.add_argument('--var_fn', type=str, default=None,
                        help=SUPPRESS)
    args = parser.parse_args()

    CreateTensorFullAlignment(args)


if __name__ == "__main__":
    main()