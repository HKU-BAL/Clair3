import collections
import heapq
import itertools
import json
import shlex
import signal
import sys
import os

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict
from subprocess import PIPE, Popen
from shared.command_options import (
    CommandOption,
    CommandOptionWithNoValue,
    ExecuteCommand,
    command_string_from,
    command_option_from
)
from shared.utils import file_path_from, executable_command_string_from, subprocess_popen, str2bool, log_warning

from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, region_from, reference_sequence_from
import shared.param_p as param

class InstancesClass(object):
    def __init__(self):
        self.create_tensor = None

    def poll(self):
        self.create_tensor.poll()


c = InstancesClass()

def check_return_code(signum, frame):
    c.poll()
    if c.create_tensor.returncode != None and c.create_tensor.returncode != 0:
        c.compress_tensor.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if c.create_tensor.returncode == None:
        signal.alarm(5)


extended_window_size = 200
region_size = 50000
reference_region_size = region_size * 2
extend_bp = 100
reference_allele_gap = 0.9

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=sys.stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=sys.stderr, bufsize=bufsize, universal_newlines=True)

class Reference(object):
    """
    Reference region query with given reference start and reference end, we cocnat the reference base with altertive base
    with reference base to generate read-level haplotype.
    """
    def __init__(self, seq, start, reference_sequence, reference_start):
        self.start = start
        self.end = start + len(seq)
        self.seq = seq
        self.reference_sequence = reference_sequence
        self.reference_start = reference_start

    def query(self, start, end):
        return self.seq[start - self.start:end - self.start]


def all_genotypes_combination(variant, alt_dict, variant_dict):
    """
    Enumerate true variant site and candidate site genotype combination and find read. For a phased confident site, we
    only enumerate one genotype with confident flag.
    """

    if variant.type == 'candidate':
        num_ref_and_alts = len(variant.variant.alternate_bases)

        if variant.variant.start in alt_dict and alt_dict[variant.variant.start].phased_genotype:
            return 1
        elif variant.variant.start in alt_dict and alt_dict[variant.variant.start].confident_variant:
            return 2
        return (num_ref_and_alts + 1) * num_ref_and_alts / 2
    else:
        if variant.variant.start in alt_dict and alt_dict[
            variant.variant.start].phased_genotype and variant.variant.start in variant_dict:
            return 1
        return len(variant.variant.alternate_bases)

def unique_genotypes_selection(genotype_options):
    """
    Extend true two haplotypes according to the chosen genotype and only save the haplotype set with distinct haplotype,
    if two haplotypes set match in either hap1(1) == hap2(1) or hap1(1) = hap2(1), remove the duplication to reduce
    massive comparison.
    """

    genotypes_list = []
    existed_genotypes_path = set()
    for genotype_pair in itertools.product(*genotype_options):
        g1, g2 = "", ""
        for h1, h2 in genotype_pair:
            g1 += str(h1)
            g2 += str(h2)
        genotype_tupe = (g1, g2)
        genotype_tupe_reverse = (g2, g1)
        if genotype_tupe in existed_genotypes_path or genotype_tupe_reverse in existed_genotypes_path:
            continue
        existed_genotypes_path.add(genotype_tupe)
        existed_genotypes_path.add(genotype_tupe_reverse)
        genotypes_list.append(genotype_pair)
    return genotypes_list

def find_read_support(variants, ref, variant_type, max_calculate_count, variant_dict=None, read_name_info_dict=None, truths=None, alt_dict=None, no_match_found=False):
    """
    Find read-level support for each matched haplotype, we only extended the reference sequence with the alternative base,
    and discard low allele frequency systematic error.
    """
    read_seqs_counter = None
    if variant_type == 'candidate':
        all_read_set = set()
        all_read_seq_dict = defaultdict(str)
        for v in variants:
            pos = v.start
            read_name_set = alt_dict[pos].read_name_set
            all_read_set = all_read_set.union(read_name_set)
        for read_name in all_read_set:
            ref_start = ref.start
            ref_end = ref.end
            read_seq = read_name_info_dict[read_name].seq
            if not len(read_seq):
                continue
            ref_offset, alt_offset, pre_end = 0, 0, 0
            for start, end, seq in read_seq:
                if end < ref_start or start > ref_end:
                    continue
                start_offset = start - ref_start
                end_offset = end - ref_start
                if start_offset >= pre_end:
                    all_read_seq_dict[read_name] += ref.seq[pre_end:start_offset] + seq
                    pre_end = end_offset
            if pre_end < len(ref.seq):
                all_read_seq_dict[read_name] += ref.seq[pre_end:]
        read_seqs_counter = Counter(all_read_seq_dict.values())

    def extend_genotype(variants_and_genotypes, next_pos_list):

        """
        We give two iterator of two haplotype start position and update them separately, which allow haplotype extension
        more flexible when one of the haplotype has insertion or deletion, if the start position reach our provided region,
        then the extension will stop.
        """
        if next_pos_list is None or None in next_pos_list:
            pass
        if not variants_and_genotypes:
            hap1_last_pos, hap2_last_pos = next_pos_list

            rest_seq_1 = ref.query(hap1_last_pos,
                                             ref.end) if hap1_last_pos != ref.end and hap1_last_pos else ""
            rest_seq_2 = ref.query(hap2_last_pos,
                                             ref.end) if hap2_last_pos != ref.end and hap2_last_pos else ""
            yield [rest_seq_1, rest_seq_2]  # add last padding ref base
        else:
            current_variant, remaining_variants = [variants_and_genotypes[0]], variants_and_genotypes[1:]
            prefix_seqs_list, next_pos_list = find_seqs(
                current_variant, next_pos_list, ref)
            prefix_seqs_list = list(prefix_seqs_list)

            if not prefix_seqs_list or next_pos_list is None:
                pass

            for seqs in extend_genotype(remaining_variants, next_pos_list):
                yield [prefix_seqs_list[0] + seqs[0], prefix_seqs_list[1] + seqs[1]]

    def extend_genotypes(variants_and_genotypes, next_pos):

        try:
            for r in extend_genotype(variants_and_genotypes, next_pos):
                yield r
        except:
            pass

    genotypes_combinations = genotypes_combination(variants, variant_type, variant_dict, max_calculate_count, truths, alt_dict=alt_dict, no_match=no_match_found)
    genotypes_seqs_dict = defaultdict(list)
    genotypes_list = unique_genotypes_selection(genotypes_combinations)
    GT = collections.namedtuple('GT', ['variant', 'genotypes'])
    for genotypes in genotypes_list:
        variants_and_genotypes = [GT(v, g) for v, g in zip(variants, genotypes)]
        for seqs in extend_genotypes(variants_and_genotypes, [ref.start, ref.start]):
            genotypes_seqs_dict[frozenset(seqs)].append(genotypes)
    return genotypes_seqs_dict, read_seqs_counter

def remove_common_suffix(ref_base, alt_base):
    """
    For each haploid match, we simplify the reference base and alternative base and remove their common suffix characters.
    """

    min_length = min(len(ref_base) - 1, min([len(item) - 1 for item in alt_base]))  # keep at least one base
    prefix = ref_base[::-1]
    for string in alt_base:
        string = string[::-1]
        while string[:len(prefix)] != prefix and prefix:
            prefix = prefix[:len(prefix) - 1]
        if not prefix:
            break
    res_length = len(prefix)
    if res_length > min_length:
        return ref_base, alt_base
    return ref_base[:len(ref_base) - res_length], [item[:len(item) - res_length] for item in alt_base]

    return ref_base[-min_length], [item[-min_length] for item in alt_base]


def has_multi_in_truths(truths, start=None):
    for t in truths:
            if len(t.alternate_bases) > 1:
                return True
    return False

def count_combination(genotypes_combinations):
    """
    Calculate the Cartesian product required for all genotype combinations
    """
    product = 1
    for gc in genotypes_combinations:
        product *= len(gc)
    return product

def genotypes_combination(variants, variant_type, variant_dict, max_calculate_count, truths=None, alt_dict=None, no_match=False, simplfy_combination=False):

    """
    Calculate genotype combination for haplotype set generation. For a locked confident variant or candidate site, we directly
    extend with its phased genotype, while for other candidates, we need to assign with a missing flag to skip the candidate, or a
    pass flag keep the genotype for further extension.
    """

    if no_match:
        [{(0, 0)}] * len(variants)
    if variant_type == 'truth':
      output = []
      for v in variants:
        if v:
          gt = tuple(v.genotype)

          is_confident_pos = (variant_dict[v.start].confident_variant and variant_dict[v.start].phased_genotype) if v.start in variant_dict else False
          output.append({tuple(variant_dict[v.start].phased_genotype)} if is_confident_pos else {(0, 0), tuple(gt), tuple(list(gt)[::-1])})
        else:
          output.append({(-1, -1)})
      return output
    elif variant_type == 'candidate':
        output = []
        has_multi_in_truth = has_multi_in_truths(truths)

        for v in variants:
            is_confident_pos = v.start in alt_dict and alt_dict[v.start].phased_genotype
            pos_in_truths = is_confident_pos and v.start in variant_dict and variant_dict[v.start].phased_genotype
            if pos_in_truths:
                output.append({tuple(alt_dict[v.start].phased_genotype)})
            elif is_confident_pos:
                if simplfy_combination:
                    output.append({tuple(alt_dict[v.start].phased_genotype)})
                else:
                    output.append({(0, 0), tuple(alt_dict[v.start].phased_genotype)})
            else:
                gt_set = set()
                for idx_1 in range(len(v.alternate_bases) + 1):
                    for idx_2 in range(len(v.alternate_bases) + 1):
                        if simplfy_combination and not has_multi_in_truth:
                                if idx_1 != 0 and idx_2 != 0 and idx_1 != idx_2:
                                    continue
                        gt_set.add((idx_1, idx_2))
                output.append(gt_set)
        # in extra mass combination, need to simplfy low possible cases:
        if  count_combination(output) > max_calculate_count:
            if simplfy_combination:
                # skip
                return [{(0, 0)}] * len(variants)
            else:
                return genotypes_combination(variants, variant_type, variant_dict, max_calculate_count, truths, alt_dict, no_match, simplfy_combination=True)
        return output


def find_seqs(variants_and_genotypes, last_pos_list, ref):
    seqs = ["", ""]
    genotypes = [vg.genotypes for vg in variants_and_genotypes]
    variants = [vg.variant for vg in variants_and_genotypes]
    all_genotypes = [tuple([item[i] for item in genotypes]) for i in [0,1]]
    next_last_pos_list = [0,0]
    for idx, phased_genotype in enumerate(all_genotypes):
        current_seq, hap_end = build_seq(variants, phased_genotype, ref, last_pos_list[idx], None) # if overlap, merge multiple variants together
        next_last_pos_list[idx] = hap_end
        if current_seq:
            seqs[idx] = current_seq
    return seqs, next_last_pos_list


def build_seq(variants, phased_genotype, ref, pre_start, ref_end=None):

    """
    Build or extend the haplotype according to provided genotype. We marked the start position iterator of each haplotype and
    update with variant alternative base.
    """

    seqs = ""
    position = pre_start
    for variant, phased in zip(variants, phased_genotype):
        if variant.start < pre_start:
            if variant.start == pre_start - 1 and phased != 0:  # this only happen when pre pos is deletion and current pos is insertion
                ref_base = variant.reference_bases
                alt_base = variant.alternate_bases[phased - 1]
                if len(alt_base) > len(ref_base): # is an insertion
                    # print ('has insertion and deletion overlap'.format(variant.start))
                    return alt_base[1:], position
            if phased != 0: # impossible # sometimes happen in true vcf
                return None, None
            else:
                return "", pre_start # do not do anything if 0 allele
        else:
            seqs += ref.query(pre_start, variant.start)

        allele = variant.reference_bases if phased == 0 else variant.alternate_bases[phased - 1]
        if phased == 0:
            allele = allele[0]
            position = variant.start + 1
            seqs += allele # only add one ref base
        else:
            ref_base = variant.reference_bases
            alt_base = variant.alternate_bases[phased-1]
            ref_base, alt_base = remove_common_suffix(ref_base, [alt_base])
            end = variant.start + len(ref_base)
            position = end
            seqs += alt_base[0]

    return seqs, position

class ReadMatch(object):
    def __init__(self, sample_ctg_info, candidates, candidate_genotypes, truths,match_seq,truth_genotypes,match_reads_count=0):

        self.sample_ctg_info = sample_ctg_info
        self.candidates = candidates
        self.truths = truths
        self.candidate_genotypes = candidate_genotypes
        self.truth_genotypes = truth_genotypes
        self.match_reads_count = match_reads_count
        self.match_seq = sorted(match_seq)
        self.truths_pos_list = [t.start for t in truths]
        self.candidates_pos_list = [c.start for c in candidates]
        self.raw_genotypes = [tuple(v.genotype) for v in self.truths]
        self.non_variants = [int(sum(cg) == 0) for cg in self.candidate_genotypes]
        self.miss_variants_count = sum([1 if sum(gt) < sum(raw_gt) else 0 for raw_gt, gt in zip(self.raw_genotypes,self.truth_genotypes)])
        self.match_variants_count = sum([1 if item < 1 else 0 for item in self.non_variants])
        self.non_variants_count = sum(self.non_variants)
        self.match_order = (self.match_reads_count, self.miss_variants_count, self.non_variants_count, self.match_variants_count)

    def match_info(self):
        can_info, truth_info = "", ""
        for can, gt in zip(self.candidates, self.candidate_genotypes):
            gt_str = '_' + '/'.join(map(str, gt)) + ' '
            can_info += str(can.start) + '-' + can.reference_bases + '->' + '-'.join(can.alternate_bases) + gt_str

        for truth, gt in zip(self.truths, self.truth_genotypes):
            gt_str = '_' + '/'.join(map(str, gt)) + ' '
            truth_info += str(truth.start) + '-' + truth.reference_bases + '->' + '-'.join(truth.alternate_bases) + gt_str

        extro_info = ""
        if self.match_reads_count >=-6 and self.miss_variants_count > 0 and self.match_variants_count > 0:
            extro_info = '\nthis match has few read support'
        return ('ctg_info={}, read_support,miss_variants,non_variants,match_variants={}, candidate={}, truth={} {}').format(self.sample_ctg_info,
            self.match_order,
            can_info, truth_info, extro_info)

class Position(object):
    def __init__(self, pos, genotype1, genotype2, ref_base=None, alt_base=None, candidate=False, cigar_count=None,
                 confident_variant=False, depth=None, alt_list=None, af_list=None, alt_type_mapping_dict=None):
        self.pos = pos
        self.reference_bases = ref_base
        self.candidate = candidate

        if candidate == True:
            self.alternate_bases = alt_base
        else:
            self.alternate_bases = [alt_base] if ',' not in alt_base else alt_base.split(',')

        self.start = pos
        self.end = self.pos + len(ref_base)
        self.genotype = [genotype1, genotype2]
        self.cigar_count = cigar_count
        self.confident_variant = confident_variant
        self.read_name_set = set()
        self.depth = depth
        self.variant_hap_dict = defaultdict(defaultdict)
        self.phased_genotype = None
        self.hap_count_dict = defaultdict(int)
        self.alt_list = alt_list
    def update_info(self, ref_base, alt_base, genotype):
        self.reference_bases = ref_base
        self.alternate_bases = alt_base
        self.genotype = genotype


class Read(object):
    def __init__(self, hap=0):
        self.hap = hap
        self.pos_alt_dict = defaultdict(str)
        self.start = None
        self.end = None
        self.seq = []


def decode_alt_info(cigar_count, ref_base, depth, minimum_allele_gap):
    """
    Decode the input read-level alternative information
    cigar_count: each alternative base including snp, insertion and deletion of each position
    pileup_bases: pileup bases list of each read in specific candidate position from samtools mpileup 1.10
    reference_sequence: the whole reference sequence index by contig:start-end. 0-based.
    ref_base: upper reference base for cigar calculation.
    depth: depth of candidate position for calculation.
    minimum_allele_gap: default minimum allele frequency for candidate to consider as a potential true variant for unification.
    """
    alt_type_list = []  # SNP I D
    seqs = cigar_count.split(' ')
    seq_alt_bases_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
    if '*' in seq_alt_bases_dict:
        del seq_alt_bases_dict['*']
    max_del_cigar = ""
    del_list = []
    ref_represatation = ref_base
    alt_list = sorted(list(seq_alt_bases_dict.items()), key=lambda x: x[1], reverse=True)

    seq_insertion_bases_list = alt_list[:2]
    af_list = []
    for alt_type, count in seq_insertion_bases_list:
        count = int(count)
        if '*' in alt_type or '#' in alt_type:
            continue
        if count / float(depth) < minimum_allele_gap:
            continue
        af_list.append(count/ float(depth))
        if alt_type[0] == 'X':
            alt_type_list.append(alt_type[1])
        elif alt_type[0] == 'I':
            alt_type_list.append(alt_type[1:])
        elif alt_type[0] == 'D':
            if len(alt_type[1:]) > len(max_del_cigar):
                max_del_cigar = alt_type[1:]
            del_list.append(ref_base + alt_type[1:])
    new_del_list = []
    if len(max_del_cigar):
        ref_represatation = ref_base + max_del_cigar
        alt_type_list = [item + max_del_cigar for item in alt_type_list]
        for idx, item in enumerate(del_list):
            start_pos = len(item[1:])
            append_del_bases = max_del_cigar[start_pos:]
            new_del_list.append(
                ref_base + append_del_bases)  # ACG-> A, ACGTT -> A, max_del_cigar is CGTT, represent ACG-> A to ACGTT->ATT
    alt_base_list = alt_type_list + new_del_list
    return ref_represatation, alt_base_list,af_list, alt_list

def has_variant_suport(ref_base, alt_base, pos, alt_dict):
    """
    ref_base: reference base of the true varaint.
    alt_base: alternative base of the true variant
    pos: pos: candidate position for unification.
    alt_dict: dictionary (pos: pos info) which keep position level candidate reference base and alternative base information.
    return the alternative index of each candidate site if found match
    """

    alt_index = -1
    if pos not in alt_dict or not len(alt_dict[pos]):
        return alt_index
    cigar_count = alt_dict[pos]
    if len(ref_base) == 1 and len(alt_base) == 1:  # Snp
        if alt_base in cigar_count:
            alt_index = cigar_count[0][alt_base]
    elif len(ref_base) > len(alt_base):  # D
        del_base = ref_base[1:len(ref_base) - len(alt_base) + 1]
        if del_base in cigar_count:
            alt_index = cigar_count[1][alt_base]
    elif len(ref_base) < len(alt_base):  # I
        ins_base = alt_base[1:len(alt_base) - len(ref_base) + 1]
        if ins_base in cigar_count:
            alt_index = cigar_count[2][alt_base]
    return alt_index


def get_ref(ref_fn, ctg_name):
    refernce_sequences = []
    samtools_faidx_process = subprocess_popen(
        shlex.split("samtools faidx {} {}".format(ref_fn, ctg_name))
    )
    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    reference_sequence = "".join(refernce_sequences[1:])

    reference_sequence = reference_sequence.upper()
    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    reference_start = 1
    return reference_sequence, reference_start


def get_genotype(genotype):
    g1, g2 = genotype
    min_gt = min(int(g1), int(g2))
    max_gt = max(int(g1), int(g2))
    return str(min_gt) + '/' + str(max_gt)

def lock_variant(variant, truth):
    """
    Find potential locked true variant and candidate site if we consider it as a confident site, there are only exactly one
    candidate match with the true variant in a specific site, and we can further take the candidate into consideration with
    a match flag.
    """
    if truth is None:
        return None, False
    variant_ref_base = variant.reference_bases
    variant_alt_base = variant.alternate_bases
    truth_ref_base = truth.reference_bases
    truth_alt_base = truth.alternate_bases

    if len(variant_alt_base) != len(truth_alt_base):
        return None, False
    tmp_alt_list = []
    for ab in variant_alt_base:
        ref_base1, alt_base1 = remove_common_suffix(variant_ref_base, [ab])
        tmp_alt_list.append((alt_base1[0], ref_base1))
    match_index = [-1] * len(truth_alt_base)

    for t_idx, ab in enumerate(truth_alt_base):
        ref_base1, alt_base1 = remove_common_suffix(truth_ref_base, [ab])
        for a_idx, (alt_base, ref_base) in enumerate(tmp_alt_list):
            if alt_base1[0] == alt_base and ref_base1 == ref_base:
                match_index[t_idx] = a_idx
    match = sum([1 if item >= 0 else 0 for item in match_index]) == len(truth_alt_base) # can find all alt_base
    return match_index, match

def decode_variant(variant, reference_base):
    if variant == 'R':
        return 'R', 'R'
    if variant[0] == 'X':
        return reference_base, variant[1]
    elif variant[0] == 'I':
        return reference_base,  variant[1:]
    elif variant[0] == 'D':
        return reference_base + variant[1:], reference_base

def update_variant_hap_dict(alt_dict, pos, reference_sequence, reference_start, is_variant_confident, variant_dict, allele_gap, platform):
    """
    For a phased alignment, the candidates are easier to lock as confident if the signal exists strongly in one side and have confident
    match with true variant.
    """
    phased_genotype = [-1,-1]
    reference_base = reference_sequence[pos- reference_start]
    variant_hap_dict = alt_dict[pos].variant_hap_dict
    if not len(variant_hap_dict):
        return None
    hap_count_dict = alt_dict[pos].hap_count_dict
    variant_ref_base = alt_dict[pos].reference_bases
    variant_alt_base = alt_dict[pos].alternate_bases
    for variant, hap_dict in variant_hap_dict.items():
        if variant in '*#':
            continue
        hap_0 = hap_dict[0] if 0 in hap_dict else 0
        # for illumina unification, phased information contributes less, we safely denote with reference allele gap
        if platform == 'ilmn':
            hap_total_af = (hap_0) / float(sum(list(hap_count_dict.values())))
            if variant not in 'R*' and hap_total_af > reference_allele_gap and is_variant_confident and variant_dict[pos].genotype == [1, 1]:
                phased_genotype = [1, 1]
                return phased_genotype
            if -1 in phased_genotype:
                return None
            return phased_genotype

        hap_1 = hap_dict[1] if 1 in hap_dict else 0
        hap_2 = hap_dict[2] if 2 in hap_dict else 0
        hap_0 = hap_0 if hap_0 > 3 else 0
        hap_1 = hap_1 if hap_1 > 3 else 0
        hap_2 = hap_2 if hap_2 > 3 else 0

        hap_total_af = (hap_0 + hap_1 + hap_2) / float(sum(list(hap_count_dict.values())))
        if variant not in 'R*#' and hap_total_af > 1 - allele_gap / 2 and is_variant_confident and variant_dict[pos].genotype== [1, 1]:
            phased_genotype = [1, 1]
            return phased_genotype
        hap_1_af = hap_1 / float(hap_count_dict[1]) if 1 in hap_count_dict else 0.0
        hap_2_af = hap_2 / float(hap_count_dict[2]) if 2 in hap_count_dict else 0.0

        if variant == 'R':
            if hap_1_af >= 1 - allele_gap:
                phased_genotype[0] = 0
            if hap_2_af >= 1 - allele_gap:
                phased_genotype[1] = 0
            continue
        ref_base, alt_base = decode_variant(variant, reference_base)

        alt_index = -1
        for ab_idx, ab in enumerate(variant_alt_base):
            ref_base1, alt_base1 = remove_common_suffix(variant_ref_base, [ab])
            if alt_base1[0] == alt_base and ref_base1 == ref_base:
                alt_index = ab_idx
                break
        if alt_index == -1:
            continue
        if hap_1_af >= 1 - allele_gap * 2:
            phased_genotype[0] = alt_index + 1
        if hap_2_af >= 1 - allele_gap * 2:
            phased_genotype[1] = alt_index + 1

    if -1 in phased_genotype:
        return None
    return phased_genotype


def match_alt_base(alt_list, ref_base, alt_base):
    if not len(alt_list) or (len(alt_list) == 1 and 'R' in alt_list):
        return False
    alt_set = set([item[0] for item in alt_list])

    for ab in alt_base:
        rb, ab = remove_common_suffix(ref_base, [ab])
        if len(rb) == len(ab[0]): #snp
            ab = 'X' + ab[0]
            if ab in alt_set:
                return True
        elif len(rb) < len(ab[0]): # insertion
            ab = 'I' + ab[0]
            if ab in alt_set:
                return True
        elif len(rb) > len(ab[0]):
            ab = 'D' + rb[1:]
            if ab in alt_set:
                return True
    return False

def check_confident_match(candidates, truths):

    """
    Double check whether the candidate site match the representation in the truth variant site in reference base,
    alternative base, genotype and start position.
    """

    if len(candidates) != len(truths):
        return False
    all_candidate_positions = set([c.start for c in candidates])
    for truth in truths:
        if truth.start not in all_candidate_positions:
            return False
        for candidate in candidates:
            if candidate.start == truth.start:
                if candidate.reference_bases != truth.reference_bases or candidate.alternate_bases != truth.alternate_bases:
                    return False
    return True

def split_variants_truths(candidates,
                          truths,
                          partition_size,
                          max_candidates_distance,
                          max_calculate_count,
                          variant_dict=None,
                          alt_dict=None):
    """
    Split all candidate sites and true variant according to the start position, for true variant site, we extend the
    at least one candidate site in both two sides to aviod missing match.
    """
    INFO = collections.namedtuple('INFO', ['start', 'type', 'variant'])
    def match_max_candidate_distance(partition, variants, new_count):
        if not partition:
            return True
        n_of_type = sum(1 for g in partition if g.type == variants.type)
        if new_count >= max_calculate_count or n_of_type >= partition_size:
            if new_count >= max_calculate_count:
                print('{} exceed max calculation count'.format(new_count))
            return False
        else:
            for g in partition:
                if variants.variant.start - g.variant.end + 1 > max_candidates_distance:
                    return False

            last_par = partition[-1].variant.end
            if variants.variant.start - last_par + 1 > extend_bp:
                return False
            return True

    truths_pos_set = set([v.start for v in truths])
    sorted_variants = list(heapq.merge(
        [INFO(v.start, 'candidate', v) for v in candidates],
        [INFO(t.start, 'truth', t) for t in truths]))

    all_partitions = []
    partition = []
    product_count = 1
    for sv_idx in range(len(sorted_variants)):
        variants = sorted_variants[sv_idx]
        new_count = product_count * all_genotypes_combination(
            variants, variant_dict, alt_dict)
        if match_max_candidate_distance(partition, variants,
                                        new_count):
            partition.append(variants)
            product_count = new_count
        else:
            if variants.start == partition[-1].start and variants.type != partition[
                -1].type:  #
                # add same truths or variants together and add at least one nearby candidate
                partition.append(variants)
                if sv_idx < len(sorted_variants) - 1 and sorted_variants[sv_idx + 1].start not in truths_pos_set and \
                        sorted_variants[sv_idx + 1].start - variants.start <= extend_bp:
                    partition.append(sorted_variants[sv_idx + 1])
                all_partitions.append(partition)
                partition = []
                product_count = 1
            else:
                all_partitions.append(partition)
                partition = [variants]
                product_count = all_genotypes_combination(variants, variant_dict, alt_dict)
    if partition:
        all_partitions.append(partition)

    split_partitions = []
    for partitions in all_partitions:
        candidate_partitions = []
        truth_partitions = []
        for p in partitions:
            if p.type == 'candidate':
                candidate_partitions.append(p.variant)
            elif p.type == 'truth':
                truth_partitions.append(p.variant)
        split_partitions.append([candidate_partitions, truth_partitions])
    return split_partitions

class RepresentationUnification(object):

    def __init__(self,
                 sample_name,
                 contig_name,
                 reference_sequence,
                 reference_start,
                 partition_size,
                 max_candidates_distance,
                 max_calculate_count,
                 subsample_ratio):

        self.sample_name = sample_name
        self.contig_name = contig_name
        self.subsample_ratio = subsample_ratio
        self.sample_ctg_info = '_'.join([sample_name, str(subsample_ratio), contig_name])
        self.partition_size = partition_size
        self.max_candidates_distance = max_candidates_distance
        self.max_calculate_count = max_calculate_count
        self.reference_sequence = reference_sequence
        self.reference_start = reference_start


    def get_reference_seq(self, candidates, true_variants, bufsize=50):
        all_variants = candidates + true_variants
        start = min(x.start for x in all_variants)
        end = max(x.end for x in all_variants)

        ref_bases = self.reference_sequence[start - self.reference_start - 1:end + bufsize - self.reference_start]
        return Reference(seq=ref_bases,
                         start=start - 1,
                         reference_sequence=self.reference_sequence,
                         reference_start=self.reference_start)

    def find_match_pairs(self, candidates, truths, ref, variant_dict, read_name_info_dict=None, alt_dict=None):
        no_match_found = len(candidates) == 0 or len(truths) == 0

        if no_match_found:
            can_info, truth_info = "", ""
            for can in candidates:
                gt = can.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                can_info += str(can.start) + '-' + can.reference_bases + '->' + '-'.join(can.alternate_bases) + gt_str

            for truth in truths:
                gt = truth.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                truth_info += str(truth.start) + '-' + truth.reference_bases + '->' + '-'.join(
                    truth.alternate_bases) + gt_str

            print ('[INFO] Missing match: ctg_info={}, read_support,miss_variants,non_variants,match_variants=(0, {}, {}, 0), candidate={}, truth={}'.format(self.sample_ctg_info,len(truths), len(candidates),can_info, truth_info))
            return None

        confident_match = check_confident_match(candidates, truths)
        if confident_match:
            can_info, truth_info = "", ""
            for can in candidates:
                gt = can.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                can_info += str(can.start) + '-' + can.reference_bases + '->' + '-'.join(can.alternate_bases) + gt_str

            for truth in truths:
                gt = truth.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                truth_info += str(truth.start) + '-' + truth.reference_bases + '->' + '-'.join(
                    truth.alternate_bases) + gt_str

            print (
                '[INFO] Found confident match: ctg_info={}, read_support,miss_variants,non_variants,match_variants=(None, 0, 0, {}), candidate={}, truth={}'.format(
                    self.sample_ctg_info, len(truths), can_info, truth_info))

            match_genotype = [tuple(v.genotype) for v in truths]

            return ReadMatch(
                sample_ctg_info=self.sample_ctg_info,
                candidates=truths,
                candidate_genotypes=match_genotype,
                truths=truths,
                truth_genotypes=match_genotype,
                match_seq=ref.seq,
                match_reads_count=100)

        truths_candidate_gentoypes = genotypes_combination(truths, 'truth', variant_dict, max_calculate_count,
                                                       truths, alt_dict=alt_dict, no_match=no_match_found)
        candidates_candidate_gentoypes = genotypes_combination(candidates, 'candidate', variant_dict,
                                                                              max_calculate_count,
                                                                              truths, alt_dict=alt_dict,
                                                                              no_match=no_match_found)

        truths_genotypes_list = unique_genotypes_selection(truths_candidate_gentoypes)
        candidates_genotypes_list = unique_genotypes_selection(candidates_candidate_gentoypes)

        print (len(truths_genotypes_list) * len(candidates_genotypes_list))
        if len(truths_genotypes_list) * len(candidates_genotypes_list) > self.max_calculate_count:
            return None

        truth_seqs, _ = find_read_support(
            variants=truths,
            ref=ref,
            variant_type='truth',
            max_calculate_count=self.max_calculate_count,
            variant_dict=variant_dict,
            truths=truths,
            read_name_info_dict=read_name_info_dict,
            alt_dict=alt_dict,
            no_match_found=no_match_found)

        variant_seqs, read_seqs_counter = find_read_support(
            variants=candidates,
            ref=ref,
            variant_type='candidate',
            max_calculate_count=self.max_calculate_count,
            variant_dict=variant_dict,
            truths=truths,
            read_name_info_dict=read_name_info_dict,
            alt_dict=alt_dict,
            no_match_found=no_match_found)

        matches = []
        for variant_seq, variant_genotypes in variant_seqs.items():
            if variant_seq not in truth_seqs:
                continue
            truth_seq = truth_seqs[variant_seq]
            for variant_genotype in variant_genotypes:
                match_reads_count = -sum([read_seqs_counter[seq] if read_seqs_counter and seq in read_seqs_counter else 0 for seq in variant_seq ]) # more match reads and negative count is better if smaller
                matches.append(ReadMatch(
                        sample_ctg_info=self.sample_ctg_info,
                        candidates=candidates,
                        candidate_genotypes=variant_genotype,
                        truths=truths,
                        truth_genotypes=truth_seq[0],
                        match_seq=variant_seq,
                        match_reads_count=match_reads_count))
        if not matches:
            return None
        else:
            best_matches = sorted(matches, key=lambda x: x.match_order)[0]
            print ('[INFO] Found match case:', best_matches.match_info())
            return best_matches


    def unify_label(self, variants, truths, region, ctg_start, ctg_end, all_pos, variant_dict,
                  rescue_dict=None, output_vcf_fn=None, test_pos=None, read_name_info_dict=None, alt_dict=None):
        split_start, split_end = region

        all_partitions = split_variants_truths(
            candidates=list(variants),
            truths=list(truths),
            partition_size=self.partition_size,
            max_candidates_distance=self.max_candidates_distance,
            max_calculate_count=self.max_calculate_count,
            variant_dict=variant_dict,
            alt_dict=alt_dict)

        for all_candidates, all_truths in all_partitions:
            ref = self.get_reference_seq(all_candidates, all_truths)
            match_pairs = self.find_match_pairs(candidates=all_candidates,
                                                truths=all_truths,
                                                ref=ref,
                                                variant_dict=variant_dict,
                                                read_name_info_dict=read_name_info_dict,
                                                alt_dict=alt_dict)

            if match_pairs is None:
                if not len(truths):
                    continue
                # double check to rescue true variants
                for truth in all_truths:
                    pos = truth.start
                # add missing low-confident tp position
                    if not (pos >= split_start and pos < split_end) or (ctg_start is not None and ctg_end is not None
                                                                        and not (pos >= ctg_start and pos < ctg_end)):
                        continue
                    if pos in alt_dict and pos in variant_dict:
                        ref_base = variant_dict[pos].reference_bases
                        alt_base = variant_dict[pos].alternate_bases
                        alt_list = alt_dict[pos].alt_list
                        if not match_alt_base(alt_list, ref_base, alt_base):
                            print('[INFO] {} {} miss and has no cigar support'.format(self.sample_ctg_info, pos))
                            continue
                        print('[INFO] {} {} miss by match, append to vcf'.format(self.sample_ctg_info, pos))
                        if pos in all_pos or pos in rescue_dict:
                            continue
                        ref_base = variant_dict[pos].reference_bases
                        variant = ','.join(variant_dict[pos].alternate_bases)
                        genotype_string = '/'.join(map(str, variant_dict[pos].genotype))
                        # For efficiency, we currently only compute reference base, altnertive base and genotype from GetTruth.py
                        rescue_dict[pos] = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                                self.contig_name,
                                pos,
                                ref_base,
                                variant,
                                10,
                                'PASS',
                                '.',
                                genotype_string,
                                10,
                                10,
                                0.5)
                    else:
                        print('[INFO] {} {} miss and no variant support'.format(self.sample_ctg_info, pos))
                continue
            for idx, (candidate, candidate_genotypes) in enumerate(
                    zip(match_pairs.candidates, match_pairs.candidate_genotypes)):
                pos = candidate.start

                have_miss_variants = True if sum([1 for gt in match_pairs.truth_genotypes if sum(gt) == 0]) else False

                # append a position into rescue queue if it was missed by the unification
                if sum(candidate_genotypes) == 0 and pos not in rescue_dict and have_miss_variants and pos not in variant_dict and pos in alt_dict and alt_dict[pos].phased_genotype:
                    genotype_string = '/'.join(map(str, alt_dict[pos].phased_genotype))
                    variant = ','.join(candidate.alternate_bases)
                    ref_base = candidate.reference_bases
                    rescue_dict[pos] = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                        self.contig_name, pos, ref_base, variant, 10, 'PASS', '.', genotype_string, 10, 10, 0.5)
                    continue
                if sum(candidate_genotypes) == 0:
                    continue
                if not len(candidate.alternate_bases):
                    continue
                g1, g2 = candidate_genotypes
                variant = set()
                ref_base = candidate.reference_bases
                if g1 != 0:
                    variant.add(candidate.alternate_bases[g1 - 1])
                if g2 != 0:
                    variant.add(candidate.alternate_bases[g2 - 1])
                if g1 == 0 or g2 == 0:
                    genotype_string = '0/1'
                elif g1 == g2:
                    genotype_string = '1/1'
                elif g1 != g2:
                    genotype_string = '1/2'
                ref_base, variant = remove_common_suffix(ref_base, list(variant))
                variant = ','.join(variant)
                if candidate.start in all_pos:
                    continue
                all_pos.add(pos)

                if output_vcf_fn is not None:
                    # For efficiency, we only compute reference base, altnertive base and genotype for GetTruth.py currently
                    print("%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                        self.contig_name, candidate.start, ref_base, variant, 10, 'PASS', '.', genotype_string, 10, 10, 0.5), file=output_vcf_fn)
                    if pos in rescue_dict:
                        del rescue_dict[pos]
            for idx, (pos, raw_genotype, truth_genotype) in enumerate(
                    zip(match_pairs.truths_pos_list, match_pairs.raw_genotypes,
                        match_pairs.truth_genotypes)):
                if not (pos >= split_start and pos < split_end) or (ctg_start is not None and ctg_end is not None
                                                                    and not (pos >= ctg_start and pos < ctg_end)):
                    continue

                if truth_genotype == (0, 0) and sum(raw_genotype) > 0:# miss genoytpe
                    if pos in alt_dict and pos in variant_dict:
                        ref_base = variant_dict[pos].reference_bases
                        alt_base = variant_dict[pos].alternate_bases
                        alt_list = alt_dict[pos].alt_list
                        if not match_alt_base(alt_list, ref_base, alt_base):
                            print('{} {} miss and has no cigar support'.format(self.sample_ctg_info, pos))
                            continue
                        print('{} {} miss by match, append to vcf'.format(self.sample_ctg_info, pos))
                        if pos in all_pos:
                            continue
                        all_pos.add(pos)

                        ref_base = variant_dict[pos].reference_bases
                        variant = ','.join(variant_dict[pos].alternate_bases)
                        genotype_string = '/'.join(map(str, variant_dict[pos].genotype))

                        if output_vcf_fn is not None:
                            rescue_dict[pos] = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                                self.contig_name, pos, ref_base, variant, 10, 'PASS', '.', genotype_string, 10, 10, 0.5)


def parse_candidates_file(candidate_details_fn, contig_name=None):
    candidate_details_list = []
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % candidate_details_fn))
    for row in unzip_process.stdout:
        candidate_details_list.append(row)
    return candidate_details_list


def parse_line(row):
    row = row.strip().split('\t')  # ['chr_pos', 'depth', 'cigar_count']
    chr_pos, depth, var_read_json = row[:3]
    ctg_name, pos = chr_pos.split()
    pos, depth = int(pos), int(depth)
    return (ctg_name, pos, depth, var_read_json)


def UnifyRepresentation(args):

    """
    Representation Unification algorithm main function, this algorithm aims to unify variant representation
    between training material and true variant set.
    All candidate sites with sufficient read support and over a certain allele frequency were selected, of which the
    same variant information with true set was locked as confident candidate sites. Secondly, for each remaining
    candidate site and true variant site, a matching or missing flag was assigned. We build  haplotype pairs based on
    all possible flag combinations. For each fully matching pair, we use the match pair with the most support reads as
    the final best match pair. For the remaining unmatched sites, we decrease allele frequency to further seek remaining
    candidate sites. We rescue the candidate sites matching true variant type. In the end, the unified VCF consists of
    locked variant sites, unified variant sites, and rescued variant sites.
    """
    sample_name = args.sampleName
    var_fn = args.var_fn  # true vcf var
    candidate_details_fn = args.candidate_details_fn
    contig_name = args.ctgName
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    bed_fn = args.bed_fn
    is_confident_bed_file_given = bed_fn is not None
    partition_size = args.partition_size
    minimum_allele_gap = args.minimum_allele_gap
    max_candidates_distance = args.max_candidates_distance
    global max_calculate_count
    max_calculate_count = args.max_calculate_count
    subsample_ratio = args.subsample_ratio
    platform = args.platform
    chunk_id = args.chunk_id
    chunk_num = args.chunk_num

    global test_pos
    test_pos = None

    alt_dict = defaultdict()
    read_name_info_dict = defaultdict(Read)

    if candidate_details_fn is None:
        basedir = os.path.dirname(__file__)
        CTFA_Bin = basedir + "/../clair3.py CreateTensorFullAlignment"
        pypyBin = executable_command_string_from(args.pypy, exit_on_not_found=True)
        bam_fn = file_path_from(args.bam_fn, exit_on_not_found=True)
        ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
        vcf_fn = file_path_from(args.vcf_fn)
        extend_bed = file_path_from(args.extend_bed)
        min_af = args.min_af
        ctgStart, ctgEnd = None, None
        if ctg_start is not None and ctg_end is not None and int(ctg_start) <= int(ctg_end):
            ctgStart = CommandOption('ctgStart', ctg_start)
            ctgEnd = CommandOption('ctgEnd', ctg_end)
        chunkId, chunkNum = None, None
        if chunk_id is not None and chunk_num is not None and int(chunk_id) <= int(chunk_num):
            chunkId = CommandOption('chunk_id', chunk_id)
            chunkNum = CommandOption('chunk_num', chunk_num)

        create_tensor_command_options = [
            pypyBin,
            CTFA_Bin,
            CommandOption('bam_fn', bam_fn),
            CommandOption('ref_fn', ref_fn),
            CommandOption('vcf_fn', vcf_fn),
            CommandOption('ctgName', contig_name),
            CommandOption('platform', platform),
            CommandOption('bed_fn', bed_fn),
            CommandOption('extend_bed', extend_bed),
            ctgStart,
            ctgEnd,
            chunkId,
            chunkNum,
            CommandOptionWithNoValue('unify_repre'),
            CommandOptionWithNoValue('phasing_info_in_bam'),
            CommandOption('unify_repre_fn', 'PIPE')
        ]
        if min_af is not None:
            create_tensor_command_options.append(CommandOption('min_af', min_af))
    else:
        candidate_details_list = []
        if os.path.exists(candidate_details_fn):
            candidate_details_list = parse_candidates_file(candidate_details_fn, contig_name)
        else:
            directory, prefix = os.path.split(candidate_details_fn)
            for f in os.listdir(directory):
                if not f.startswith(prefix):
                    continue
                candidate_details_list += parse_candidates_file(os.path.join(directory, f), contig_name)

        candidate_details_list = sorted(candidate_details_list, key=lambda x: x[0])

        if not len(candidate_details_list):
            return

        if ctg_start is None or ctg_end is None:
            alt_start = parse_line(candidate_details_list[0])[1]
            alt_end = parse_line(candidate_details_list[-1])[1]
            chunk_id = args.chunk_id - 1  # 1-base to 0-base
            chunk_num = args.chunk_num
            chunk_size = (alt_end - alt_start) // chunk_num + 1
            ctg_start = alt_start + chunk_size * chunk_id
            ctg_end = ctg_start + chunk_size

    is_ctg_name_given = contig_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    ref_regions = []
    reference_start = 1
    if is_ctg_range_given:
        reference_start, reference_end = ctg_start - reference_region_size, ctg_end + reference_region_size
        reference_start = 1 if reference_start < 1 else reference_start
        ref_regions.append(region_from(ctg_name=contig_name, ctg_start=reference_start, ctg_end=reference_end))
    elif is_ctg_name_given:
        ref_regions.append(region_from(ctg_name=contig_name))
        reference_start = 1

    reference_sequence = reference_sequence_from(
        samtools_execute_command='samtools',
        fasta_file_path=args.ref_fn,
        regions=ref_regions
    )

    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (var_fn)))
    variant_dict = defaultdict(list)
    for row in unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        if ctg_start is not None and ctg_end is not None and \
            (pos < ctg_start - extended_window_size or pos > ctg_end + extended_window_size):
            continue
        ref_base = columns[2]
        alt_base = columns[3]
        genotype1 = int(columns[4])
        genotype2 = int(columns[5])
        variant_dict[pos] = Position(pos=pos,
                                     ref_base=ref_base,
                                     alt_base=alt_base,
                                     genotype1=genotype1,
                                     genotype2=genotype2)

    if candidate_details_fn is None:
        try:
            c.create_tensor = subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options))
            )
            candidate_source = c.create_tensor.stdout

            signal.signal(signal.SIGALRM, check_return_code)
            signal.alarm(2)
        except Exception as e:
            print(e, file=sys.stderr)
            sys.exit("Failed to start required processes. Exiting...")
    else:
        candidate_source = candidate_details_list

    for row in candidate_source:
        ctg_name, pos, depth, var_read_json = parse_line(row)
        if contig_name != ctg_name:
            continue
        if ctg_start is not None and ctg_end is not None and \
            (pos < ctg_start - extended_window_size or pos > ctg_end + extended_window_size):
            continue

        var_read_dict = json.loads(var_read_json)
        if not len(var_read_dict):
            continue

        cigar_count = ' '.join([' '.join([item, str(len(var_read_dict[item].split(' ')))]) for item in var_read_dict.keys()])
        ref_base = reference_sequence[pos - reference_start]
        pos_in_truths = pos in variant_dict
        ref_base, alt_base, af_list,alt_list = decode_alt_info(cigar_count=cigar_count,
                                                               ref_base=ref_base,
                                                               depth=depth,
                                                               minimum_allele_gap=minimum_allele_gap)

        alt_dict[pos] = Position(pos=pos,
                                 ref_base=ref_base,
                                 alt_base=alt_base,
                                 genotype1=-1,
                                 genotype2=-1,
                                 candidate=True,
                                 depth=depth,
                                 alt_list=alt_list)

        for variant, read_str in var_read_dict.items():
            read_list = read_str.split(' ')
            for read_name in read_list:
                read_name, hap = read_name[:-2], read_name[-1]
                if read_name not in read_name_info_dict or read_name_info_dict[read_name].hap == 0 and hap != 0:
                    read_name_info_dict[read_name].hap = int(hap)

                read_hap = read_name_info_dict[read_name].hap if read_name in read_name_info_dict else 0
                if read_hap in alt_dict[pos].variant_hap_dict[variant]:
                    alt_dict[pos].variant_hap_dict[variant][read_hap] += 1
                else:
                    alt_dict[pos].variant_hap_dict[variant][read_hap] = 1
                alt_dict[pos].hap_count_dict[read_hap] += 1
                alt_dict[pos].read_name_set.add(read_name)
                read_name_info_dict[read_name].pos_alt_dict[pos] = variant

        match_index, is_variant_confident = lock_variant(alt_dict[pos], variant_dict[pos] if pos_in_truths else None)
        if is_variant_confident:
            variant_dict[pos].confident_variant = match_index
        alt_dict[pos].phased_genotype = update_variant_hap_dict(alt_dict=alt_dict,
                                                                pos=pos,
                                                                reference_sequence=reference_sequence,
                                                                reference_start=reference_start,
                                                                is_variant_confident=is_variant_confident,
                                                                variant_dict=variant_dict,
                                                                allele_gap=minimum_allele_gap,
                                                                platform=platform)
        # lock the candidate if it has meet the phased_genotype requirement and have a exactly one match true variant site
        if alt_dict[pos].phased_genotype and pos_in_truths and is_variant_confident:
            if alt_dict[pos].phased_genotype.count(0) != variant_dict[pos].genotype.count(0) or (sum(variant_dict[pos].genotype) == 3 and sum(alt_dict[pos].phased_genotype) != 3):
                # skip wrong genotype
                alt_dict[pos].phased_genotype = None
            variant_dict[pos].reference_bases = alt_dict[pos].reference_bases
            variant_dict[pos].alternate_bases = alt_dict[pos].alternate_bases
            variant_dict[pos].phased_genotype = alt_dict[pos].phased_genotype


    if is_confident_bed_file_given:
        tree = bed_tree_from(bed_fn, contig_name=contig_name)
    for read_name, read in read_name_info_dict.items():
        if not len(read_name_info_dict[read_name].pos_alt_dict):
            continue
        for pos, alt_base in read_name_info_dict[read_name].pos_alt_dict.items():
            read.start = min(read.start, pos) if read.start is not None else pos
            if alt_base[0] == 'X':
                read.seq.append((pos, pos+1, alt_base[1]))
            elif alt_base[0] == 'I':
                read.seq.append((pos, pos+1, alt_base[1:]))
            elif alt_base[0] == 'D':
                del_length = len(alt_base[1:])
                read.seq.append((pos, pos+del_length+1, reference_sequence[pos - reference_start]))
            else: #"R"
              continue
        read.end = max([item[1] for item in read.seq]) if len(read.seq) else None

    if candidate_details_fn is None:
        c.create_tensor.stdout.close()
        c.create_tensor.wait()
        signal.alarm(0)

    if not len(alt_dict) or not len(variant_dict):
        return

    region_start = ctg_start if ctg_start else min(alt_dict.keys())
    region_end = ctg_end if ctg_end else max(alt_dict.keys())
    rescue_dict = defaultdict()
    all_pos = set()

    output_vcf_fn = None
    if args.output_vcf_fn:
      output_vcf_fn = open(args.output_vcf_fn, "w")

      def output(string_value):
          print(string_value, file=output_vcf_fn)

      from textwrap import dedent

      output(dedent("""\
            ##fileformat=VCFv4.2
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
            ##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
            ##ALT=<ID=DEL,Description="Deletion">
            ##ALT=<ID=INS,Description="Insertion of novel sequence">
            ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
            ##INFO=<ID=P,Number=0,Type=Flag,Description="Whether calling result from pileup calling">
            ##INFO=<ID=F,Number=0,Type=Flag,Description="Whether calling result from full alignment calling">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
            ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1)">"""
                    ))
      output('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name))

    # For pacbio hifi platform, we select larger candidate distance for better performance
    if platform == 'hifi':
        max_candidates_distance = 200
        partition_size = 20

    for split_idx in range((region_end - region_start) // region_size + 1):
        split_start = region_start + split_idx * region_size
        split_end = split_start + region_size
        extend_split_start = split_start - extend_bp
        extend_split_end = split_end + extend_bp

        #bed region include last pos
        variants = sorted([(item, alt_dict[item]) for item in alt_dict.keys() if
                           item >= extend_split_start and item < extend_split_end and len(alt_dict[item].alternate_bases) and is_region_in(
                               tree=tree,
                               contig_name=contig_name,
                               region_start=item-2,
                               region_end=alt_dict[item].end + 2)], key=lambda x: x[0])
        variants = [item[1] for item in variants]

        truths = sorted([(item, variant_dict[item]) for item in variant_dict.keys() if
                         item >= extend_split_start and item < extend_split_end and is_region_in(
                             tree=tree,
                             contig_name=contig_name,
                             region_start=item-2,
                             region_end=variant_dict[item].end + 2)], key=lambda x: x[0])
        truths = [item[1] for item in truths]

        if not len(variants) and not len(truths):
            continue
        RU = RepresentationUnification(
            sample_name=sample_name,
            contig_name=contig_name,
            reference_sequence=reference_sequence,
            reference_start=reference_start,
            partition_size=partition_size,
            max_candidates_distance=max_candidates_distance,
            max_calculate_count=max_calculate_count,
            subsample_ratio=subsample_ratio)

        RU.unify_label(variants=variants,
                         truths=truths,
                         region=(split_start, split_end),
                         ctg_start=ctg_start,
                         ctg_end=ctg_end,
                         all_pos=all_pos,
                         variant_dict=variant_dict,
                         rescue_dict=rescue_dict,
                         output_vcf_fn=output_vcf_fn,
                         test_pos=test_pos,
                         read_name_info_dict=read_name_info_dict,
                         alt_dict=alt_dict)

    if not len(rescue_dict):
        return
    if output_vcf_fn is not None:
        for pos, vcf_info in rescue_dict.items():
            print(vcf_info, file=output_vcf_fn)
        output_vcf_fn.close()

    if os.path.exists(args.output_vcf_fn):
        for row in open(args.output_vcf_fn, 'r'):
            if row[0] != '#':
                return
        os.remove(args.output_vcf_fn)
        print("[INFO] No vcf output for file {}, remove empty file".format(args.output_vcf_fn))

def main():
    parser = ArgumentParser(description="Representation unification for candidate site and true variant")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--var_fn', type=str, default=None,
                        help="Truth variants list input from GetTruth.py")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--candidate_details_fn', type=str, default=None,
                        help="Read-level candidate details file, default: %(default)s")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="VCF output filename or stdout if not set,default: %(default)s")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    # options for advanced users
    parser.add_argument('--max_candidates_distance', type=int, default=100,
                        help="EXPERIMENTAL: Maximum distance between subsequent variants within a group")

    parser.add_argument('--max_calculate_count', type=int, default=10000,
                        help="EXPERIMENTAL: Maximum calculation times for chunk window ")

    parser.add_argument('--partition_size', type=int, default=15,
                        help="EXPERIMENTAL: Maximum variants in per group size")

    parser.add_argument('--minimum_allele_gap', type=int, default=0.15,
                        help="EXPERIMENTAL: Minimum allele gap filtering candidate path generation")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, will choose candidate +/- 1 or +/- 2. Use together with gen4Training. default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, will choose candidate +/- 1 or +/- 2. Use together with gen4Training. default: %(default)s")

    parser.add_argument('--extend_bed', type=str, default=None,
                        help=SUPPRESS)

    # options for internal process control
    ## Subsample ratio tag for sub-sampled BAM file
    parser.add_argument('--subsample_ratio', type=int, default=1000,
                        help=SUPPRESS)

    ## Test in specific candidate site. Only use for analysis
    parser.add_argument('--test_pos', type=int, default=0,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--min_af', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--bam_fn', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--pypy', type=str, default="pypy3",
                        help=SUPPRESS)

    args = parser.parse_args()
    UnifyRepresentation(args)


if __name__ == "__main__":
    main()
