import sys
import os
import shlex
import ctypes
import re
from subprocess import PIPE
from os.path import isfile
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import shared.param_f as param
from shared.utils import file_path_from, subprocess_popen, reference_sequence_from, \
    IUPAC_base_to_ACGT_base_dict as BASE2ACGT, IUPAC_base_to_num_dict as BASE2NUM

from shared.interval_tree import bed_tree_from
from shared.intervaltree.intervaltree import IntervalTree

min_dbg_mapping_quality = min_dbg_base_quality = 20
region_expansion_in_bp = expand_align_ref_region = 20
min_windows_distance = expand_align_ref_region * 4
max_window_size = max_region_reads_num = 1000

# using 5 charaters for store long read name
CHAR_STR = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$%&()*+./:;<=>?[]^`{|}~"
L_CHAR_STR = len(CHAR_STR)
EXP = 5
T_READ_NAME = L_CHAR_STR ** EXP
L_CHAR_STR_EXP = [L_CHAR_STR ** i for i in range(EXP - 1, 0, -1)]

realigner_mod = os.path.join(*(os.path.split(__file__)[:-1] + ('realign/realigner',)))
dbg_mod = os.path.join(*(os.path.split(__file__)[:-1] + ('realign/debruijn_graph',)))

realigner = ctypes.cdll.LoadLibrary(realigner_mod)
dbg = ctypes.cdll.LoadLibrary(dbg_mod)


class StructPointer(ctypes.Structure):
    _fields_ = [("position", ctypes.c_int * max_region_reads_num),
                ("cigar_string", ctypes.c_char_p * max_region_reads_num),
                ]


class DBGPointer(ctypes.Structure):
    _fields_ = [("consensus_size", ctypes.c_int),
                ("consensus", ctypes.c_char_p * 200),
                ]


#Read class for storing read information
cigar_indel_re = r"(\d+)(D)"
cigarRe = r"(\d+)([MIDNSHP=X])"
graph_min_mapping_quality = 14
def get_len(seq, cigar):
    if 'D' not in cigar:
        return len(seq)
    indel_length = 0
    for m in re.finditer(cigar_indel_re, cigar):
        indel_length += int(m.group(1))
    return len(seq) + indel_length


class Read(object):
    def __init__(self, read_start, seq, cigar, mapping_quality, base_quality, strand, raw_base_quality=None,
                 unalign=False, read_name=None, read_id=None, flag=None, RNEXT=0, PNEXT=0, TLEN=0, phasing=None):
        self.read_start = read_start
        self.cigar = cigar
        self.mapping_quality = mapping_quality
        self.seq = seq
        self.base_quality = base_quality
        self.read_id = read_id
        self.read_end = self.read_start + get_len(seq, cigar)
        self.strand = strand
        self.graph_mq = True if self.mapping_quality >= graph_min_mapping_quality else False
        self.raw_base_quality = raw_base_quality
        self.read_name = read_name
        self.region = {}
        self.region_cigar = None
        self.region_start = None
        self.flag = str(flag)
        self.RNEXT = RNEXT
        self.PNEXT = PNEXT
        self.TLEN = PNEXT
        self.test_pos = None
        self.best_cigar = cigar
        self.best_pos = read_start
        self.best_align_score = None
        self.phasing = phasing

    def set_realign_flag(self):
        self.unalign = True

    def count_align_score(self, cigar):
        score = 0
        for m in re.finditer(cigarRe, cigar):
            l, op, = int(m.group(1)), m.group(2)
            if op in 'MX=S':
                continue
            elif op in 'ID':
                score += l
        return score

    def set_realignment_info(self, region_start, realignment_cigar, realignment_start):
        realignment_cigar = realignment_cigar.replace('X', 'M')
        if realignment_cigar == self.cigar and realignment_start == self.read_start:
            return

        if self.best_align_score and realignment_cigar == self.best_cigar and realignment_start == self.best_pos:
            return
        realignment_align_score = self.count_align_score(realignment_cigar)
        if not self.best_align_score or realignment_align_score >= self.best_align_score:
            self.best_cigar = realignment_cigar
            self.best_pos = realignment_start
            self.best_align_score = realignment_align_score

    def decode_region(self, region_str):
        if region_str == '-' or '-' not in region_str:
            return
        region_str = region_str.rstrip().split('_')
        for region in region_str:
            region, cigar, pos = region.split('-')
            region, pos = int(region), int(pos)
            self.region[region] = [cigar, pos]



def byte(x):
    return bytes(x, encoding="utf8")


def find_max_overlap_index(query_region, search_regions):
    def overlap_length(region1, region2):
        return max(0, (min(region1[1], region2[1]) - max(region1[0], region2[0])))

    overlap_lengths = [overlap_length(query_region, search_region) for search_region in search_regions]
    argmax = max(range(len(search_regions)), key=lambda idx: overlap_lengths[idx])
    return None if overlap_lengths[argmax] == 0 else argmax


def get_reference_seq(sequence, start, end, reference_start_0_based):
    if end < start:
        end, start = start, end
    return sequence[start - reference_start_0_based: end - reference_start_0_based]


def phredscore2raw_score(qual):
    return ord(qual) - 33


def evc_base_from(base):
    return base if base == "N" else BASE2ACGT[base]


def region_from(ctg_name, ctg_start=None, ctg_end=None):
    """
    1-based region string [start, end]
    """
    if ctg_name is None:
        return ""
    if (ctg_start is None) != (ctg_end is None):
        return ""

    if ctg_start is None and ctg_end is None:
        return "{}".format(ctg_name)
    return "{}:{}-{}".format(ctg_name, ctg_start, ctg_end)


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def get_halpotype_tag(samtools_view_columns):
    found_hp_tag = False
    tag = [c for c in samtools_view_columns if 'HP:i:' in c]
    if not len(tag) or len(tag[0]) < 6 or not tag[0][5].isdigit():
        return None
    return tag[0][5]


def is_too_many_soft_clipped_bases_for_a_read_from(CIGAR):
    soft_clipped_bases = 0
    total_alignment_positions = 0

    advance = 0
    for c in str(CIGAR):
        if c.isdigit():
            advance = advance * 10 + int(c)
            continue
        if c == "S":
            soft_clipped_bases += advance
        total_alignment_positions += advance
        advance = 0

    # skip a read less than 55% aligned
    return 1.0 - float(soft_clipped_bases) / (total_alignment_positions + 1) < 0.55


def samtools_view_generator_from(samtools_view_process, aligned_reads, pileup, ctg_name, reference_sequence,
                                 reference_start_0_based, header):
    CHUNK_SIZE = param.realign_chunk_size
    chunk_start, chunk_end = None, None
    rs_idx = -1
    for row_id, row in enumerate(samtools_view_process.stdout):
        if row[0] == '@':
            header.append(row)
            continue
        columns = row.strip().split()
        RNAME = columns[2]
        if RNAME != ctg_name:
            continue

        read_name = columns[0]
        FLAG = int(columns[1])
        POS = int(columns[3]) - 1  # switch from 1-base to 0-base to match sequence index
        MAPQ = int(columns[4])
        CIGAR = columns[5]
        SEQ = columns[9].upper()  # uppercase for SEQ (regexp is \*|[A-Za-z=.]+)
        RNEXT = columns[6]
        PNEXT = columns[7]
        TLEN = columns[8]
        reference_position = POS
        query_position = 0
        raw_base_quality = columns[10]
        QUAL = [phredscore2raw_score(item) for item in raw_base_quality]
        STRAND = (16 == (FLAG & 16))
        HP_TAG = get_halpotype_tag(columns[11:])
        read_name += "_" + str(int(STRAND))  # distinguish two strand
        read_name, rs_idx = simplfy_read_name(rs_idx)
        if chunk_start is None:
            chunk_start = POS
            chunk_end = chunk_start + CHUNK_SIZE
        if POS >= chunk_end + region_expansion_in_bp:
            yield chunk_start, chunk_end
            chunk_start += CHUNK_SIZE
            chunk_end += CHUNK_SIZE

        read = Read(read_start=POS,
                    seq=SEQ,
                    cigar=CIGAR,
                    mapping_quality=MAPQ,
                    base_quality=QUAL,
                    strand=STRAND,
                    raw_base_quality=raw_base_quality,
                    read_name=read_name,
                    flag=FLAG,
                    PNEXT=PNEXT,
                    RNEXT=RNEXT,
                    TLEN=TLEN,
                    phasing=HP_TAG)

        if CIGAR == "*" or is_too_many_soft_clipped_bases_for_a_read_from(CIGAR):
            continue

        aligned_reads[read_name] = read
        if MAPQ < min_dbg_mapping_quality:
            continue
        advance = 0
        for c in str(CIGAR):
            if c.isdigit():
                advance = advance * 10 + int(c)
                continue
            if c == '=':
                reference_position += advance
                query_position += advance
            elif c == "M" or c == 'X':
                for _ in range(advance):
                    if QUAL[query_position] >= min_dbg_base_quality:
                        reference_base = reference_sequence[reference_position - reference_start_0_based]  # 0 base
                        query_base = SEQ[query_position]
                        if reference_base in 'ACGT' and query_base != reference_base:
                            pileup[reference_position]['X'] += 1
                    reference_position += 1
                    query_position += 1

            elif c == "I" or c == 'S':
                pre_base = reference_sequence[reference_position - reference_start_0_based - 1]
                ins_base_quality = QUAL[query_position: query_position + advance]
                out_of_region = reference_position < chunk_start - region_expansion_in_bp or reference_position > chunk_end + region_expansion_in_bp
                if not out_of_region and pre_base in 'ACGT' and (
                        sum([True for bq in ins_base_quality if bq < min_dbg_base_quality]) == 0):
                    # skip the bad seq
                    start = reference_position - advance
                    end = reference_position + advance
                    for ins_idx in range(start, end):
                        pileup[ins_idx]["X"] += 1

                # insertion consumes query
                query_position += advance

            elif c == "D":
                out_of_region = reference_position < chunk_start - region_expansion_in_bp or reference_position > chunk_end + region_expansion_in_bp
                pre_base = reference_sequence[reference_position - reference_start_0_based - 1]  # 0-base
                if not out_of_region and pre_base in 'ACGT':
                    start = reference_position
                    end = reference_position + advance
                    for ins_idx in range(start, end):
                        pileup[ins_idx]["X"] += 1
                # deletion consumes reference
                reference_position += advance
            # reset advance
            advance = 0

    yield chunk_start, chunk_end
    yield None, None


def simplfy_read_name(rs_idx):
    rs_idx = (rs_idx + 1) % T_READ_NAME
    save_read_name = ""
    div_num = rs_idx
    for div_exp in L_CHAR_STR_EXP:
        save_read_name += CHAR_STR[div_num // div_exp]
        div_num = div_num % div_exp
    if EXP != 1:
        save_read_name += CHAR_STR[div_num % L_CHAR_STR]
    return save_read_name, rs_idx


def reads_realignment(args):
    bed_file_path = args.full_aln_regions
    extend_bed = args.extend_bed
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    samtools_execute_command = args.samtools
    bam_file_path = args.bam_fn
    minMQ = args.minMQ
    min_coverage = args.minCoverage
    is_bed_file_given = bed_file_path is not None
    is_ctg_name_given = ctg_name is not None
    read_fn = args.read_fn

    global test_pos
    test_pos = None
    if is_bed_file_given:
        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_file_path)))
        candidate_file_path_output = candidate_file_path_process.stdout

        ctg_start, ctg_end = float('inf'), 0
        for row in candidate_file_path_output:
            row = row.rstrip().split('\t')
            if row[0] != ctg_name: continue
            position = int(row[1]) + 1
            end = int(row[2]) + 1
            ctg_start = min(position, ctg_start)
            ctg_end = max(end, ctg_end)

        candidate_file_path_output.close()
        candidate_file_path_process.wait()

    if chunk_id is not None:
        fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
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

    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None

    # 1-based regions [start, end] (start and end inclusive)
    ref_regions = []
    reads_regions = []
    reference_start, reference_end = None, None

    if is_ctg_range_given:
        extend_start = ctg_start - max_window_size
        extend_end = ctg_end + max_window_size
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

    tree = bed_tree_from(bed_file_path=bed_file_path)
    if is_bed_file_given and ctg_name not in tree:
        sys.exit("[ERROR] ctg_name({}) not exists in bed file({}).".format(ctg_name, bed_file_path))

    bed_option = ' -L {}'.format(extend_bed) if extend_bed else ""
    bed_option = ' -L {}'.format(bed_file_path) if is_bed_file_given else bed_option
    mq_option = ' -q {}'.format(minMQ) if minMQ > 0 else ""
    samtools_view_command = "{} view -h {} {}".format(samtools_execute_command, bam_file_path,
                                                      " ".join(reads_regions)) + mq_option + bed_option
    samtools_view_process = subprocess_popen(
        shlex.split(samtools_view_command)
    )

    if read_fn and read_fn == 'PIPE':
        save_file_fp = TensorStdout(sys.stdout)
    elif read_fn:
        save_file_fp = subprocess_popen(shlex.split("{} view -bh - -o {}".format(samtools_execute_command, read_fn + (
            '.{}_{}'.format(ctg_start, ctg_end) if is_ctg_range_given and not test_pos else ""))), stdin=PIPE,
                                        stdout=PIPE)

    reference_start_0_based = 0 if reference_start is None else (reference_start - 1)

    header = []
    add_header = False
    aligned_reads = defaultdict()
    pileup = defaultdict(lambda: {"X": 0})
    samtools_view_generator = samtools_view_generator_from(samtools_view_process=samtools_view_process,
                                                           aligned_reads=aligned_reads,
                                                           pileup=pileup,
                                                           ctg_name=ctg_name,
                                                           reference_sequence=reference_sequence,
                                                           reference_start_0_based=reference_start_0_based,
                                                           header=header)
    pre_aligned_reads = defaultdict()

    while True:
        chunk_start, chunk_end = next(samtools_view_generator)
        if chunk_start is None:
            break
        if not add_header:
            save_file_fp.stdin.write(''.join(header))
            add_header = True

        variant_allele_list = [[position, pileup[position]["X"]] for position in list(pileup.keys())]
        candidate_position_list = [(position, support_allele_count) for position, support_allele_count in
                                   variant_allele_list if
                                   support_allele_count >= min_coverage and position >= chunk_start - region_expansion_in_bp - 1 and position <= chunk_end + region_expansion_in_bp - 1]
        candidate_position_list.sort(key=(lambda x: x[0]))

        if not len(aligned_reads) or not len(candidate_position_list):
            continue
        if len(pre_aligned_reads):  # update the read in previous chunk
            for read_name, read in pre_aligned_reads.items():
                aligned_reads[read_name] = read

        region_dict = {}
        split_region_size = max_window_size
        region_tree = IntervalTree()
        for split_idx in range((chunk_end - chunk_start) // split_region_size):
            split_start = chunk_start + split_idx * split_region_size - region_expansion_in_bp - 1
            split_end = split_start + split_region_size + region_expansion_in_bp * 2 + 1
            region_dict[(split_start, split_end)] = []
            region_tree.addi(split_start, split_end)
        for candidate_position in candidate_position_list:
            for region in region_tree.at(candidate_position[0]):
                region_dict[(region.begin, region.end)].append(candidate_position[0])

        for key, split_candidate_position_list in region_dict.items():
            start_pos, end_pos = None, None
            windows = []
            read_windows_dict = {}
            for pos in split_candidate_position_list:
                if start_pos is None:
                    start_pos = pos
                    end_pos = pos

                elif pos > end_pos + 2 * min_windows_distance:
                    temp_window = (start_pos - min_windows_distance, end_pos + min_windows_distance)
                    windows.append(temp_window)
                    read_windows_dict[temp_window] = []

                    start_pos = pos
                    end_pos = pos
                else:
                    end_pos = pos

            if start_pos is not None:
                temp_window = (start_pos - min_windows_distance, end_pos + min_windows_distance)
                windows.append(temp_window)
                read_windows_dict[temp_window] = []
            if not len(windows): continue
            windows = sorted(windows, key=lambda x: x[0])
            max_window_end = max([item[1] for item in windows])
            # #find read windows overlap_pair
            for read_name, read in aligned_reads.items():
                if read.read_start > max_window_end: continue
                argmax_window_idx = find_max_overlap_index((read.read_start, read.read_end), windows)
                if argmax_window_idx is not None:
                    read_windows_dict[windows[argmax_window_idx]].append(read_name)

            # realignment
            for window in windows:
                start_pos, end_pos = window
                if end_pos - start_pos > max_window_size:  # or (window not in need_align_windows_set):
                    continue

                ref_start = start_pos - reference_start_0_based
                ref_end = end_pos - reference_start_0_based
                ref = reference_sequence[ref_start:ref_end]
                reads = []
                low_base_quality_pos_list = []
                # pypy binding with ctypes for DBG building
                for read_name in read_windows_dict[window]:
                    read = aligned_reads[read_name]
                    if (not read.graph_mq) or read.read_start > end_pos or read.read_end < start_pos:
                        continue
                    reads.append(read.seq)
                    low_base_quality_pos_list.append(
                        ' '.join([str(bq_idx) for bq_idx, item in enumerate(read.base_quality) if int(item) < 15]))
                totoal_read_num = len(reads)
                c_ref = byte(ref)
                read_list1 = ctypes.c_char_p(byte(','.join(reads)))
                low_base_quality_pos_array = ctypes.c_char_p(byte(','.join(low_base_quality_pos_list)))

                dbg.get_consensus.restype = ctypes.POINTER(DBGPointer)
                dbg.get_consensus.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]

                dbg_p = dbg.get_consensus(ctypes.c_char_p(c_ref), read_list1, low_base_quality_pos_array,
                                          totoal_read_num)

                c_consensus, consensus_size = dbg_p.contents.consensus, dbg_p.contents.consensus_size
                consensus = [item.decode() for item in c_consensus[:consensus_size]]

                if len(consensus) == 0 or len(consensus) == 1 and consensus[0] == ref or len(
                        read_windows_dict[window]) == 0:
                    continue
                min_read_start = min([aligned_reads[item].read_start for item in read_windows_dict[window]])
                max_read_end = max([aligned_reads[item].read_end for item in read_windows_dict[window]])
                tmp_ref_start = max(0, min(min_read_start, start_pos) - expand_align_ref_region)
                tmp_ref_end = max(max_read_end, end_pos) + expand_align_ref_region

                ref_prefix = get_reference_seq(reference_sequence, tmp_ref_start, start_pos, reference_start_0_based)
                ref_center = get_reference_seq(reference_sequence, start_pos, end_pos, reference_start_0_based)
                if tmp_ref_end < end_pos:
                    continue
                ref_suffix = get_reference_seq(reference_sequence, end_pos, tmp_ref_end, reference_start_0_based)
                ref_seq = ref_prefix + ref_center + ref_suffix

                # pypy binding with ctypes for realignment
                read_name_list = []
                totoal_read_num = min(max_region_reads_num, len(read_windows_dict[window]))
                seq_list = (ctypes.c_char_p * totoal_read_num)()
                position_list = (ctypes.c_int * totoal_read_num)()
                cigars_list = (ctypes.c_char_p * totoal_read_num)()

                for read_idx, read_name in enumerate(read_windows_dict[window]):
                    read = aligned_reads[read_name]
                    if read_idx >= totoal_read_num: break
                    seq_list[read_idx] = byte(read.seq.upper())
                    position_list[read_idx] = read.read_start
                    cigars_list[read_idx] = byte(read.cigar)
                    read_name_list.append(read_name)
                haplotypes_list = [ref_prefix + cons + ref_suffix for cons in consensus]
                haplotypes = ' '.join(haplotypes_list)

                realigner.realign_reads.restype = ctypes.POINTER(StructPointer)
                realigner.realign_reads.argtypes = [ctypes.c_char_p * totoal_read_num, ctypes.c_int * totoal_read_num,
                                                    ctypes.c_char_p * totoal_read_num, ctypes.c_char_p, ctypes.c_char_p,
                                                    ctypes.c_int,
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int]

                realigner_p = realigner.realign_reads(seq_list, position_list, cigars_list,
                                                      ctypes.c_char_p(byte(ref_seq)),
                                                      ctypes.c_char_p(byte(haplotypes)), tmp_ref_start,
                                                      len(ref_prefix), len(ref_suffix), totoal_read_num)

                realign_positions, realign_cigars = realigner_p.contents.position, realigner_p.contents.cigar_string
                read_position_list = realign_positions[:totoal_read_num]
                read_cigar_list = [item.decode() for item in realign_cigars[:totoal_read_num]]

                if len(read_name_list):
                    for read_id, read_name in enumerate(read_name_list):
                        if read_cigar_list[read_id] == "" or (
                                aligned_reads[read_name].cigar == read_cigar_list[read_id] and aligned_reads[
                            read_name].read_start == read_position_list[read_id]):
                            continue
                        # update cigar and read start position
                        aligned_reads[read_name].test_pos = test_pos
                        realignment_start = read_position_list[read_id]
                        realignment_cigar = read_cigar_list[read_id].replace('X', 'M')
                        if realignment_cigar == aligned_reads[read_name].cigar and realignment_start == aligned_reads[
                            read_name].read_start:
                            continue
                        aligned_reads[read_name].set_realignment_info(split_start, read_cigar_list[read_id],
                                                                      read_position_list[read_id])

                realigner.free_memory.restype = ctypes.POINTER(ctypes.c_void_p)
                realigner.free_memory.argtypes = [ctypes.POINTER(StructPointer), ctypes.c_int]
                realigner.free_memory(realigner_p, totoal_read_num)
        # # realignment end

        if read_fn:
            sorted_key = sorted([(key, item.best_pos) for key, item in aligned_reads.items()], key=lambda x: x[1])
            for read_name, read_start in sorted_key:
                read = aligned_reads[read_name]
                if read_start < chunk_start - region_expansion_in_bp - max_window_size:  # safe distance for save reads
                    phasing_info = 'HP:i:{}'.format(read.phasing) if read.phasing else ""
                    pass
                    read_str = '\t'.join([read_name, read.flag, ctg_name, str(read_start + 1),
                                          str(read.mapping_quality), read.best_cigar, read.RNEXT, read.PNEXT, read.TLEN,
                                          read.seq,
                                          read.raw_base_quality,
                                          phasing_info])
                    save_file_fp.stdin.write(read_str + '\n')
                    del aligned_reads[read_name]
                for pile_pos in list(pileup.keys()):
                    if pile_pos < chunk_start - region_expansion_in_bp - max_window_size:
                        del pileup[pile_pos]

    if read_fn and aligned_reads:
        sorted_key = sorted([(key, item.best_pos) for key, item in aligned_reads.items()], key=lambda x: x[1])
        for read_name, read_start in sorted_key:
            read = aligned_reads[read_name]
            phasing_info = 'HP:i:{}'.format(read.phasing) if read.phasing else ""
            read_str = '\t'.join([read_name, read.flag, ctg_name, str(read_start + 1),
                                  str(read.mapping_quality), read.best_cigar, read.RNEXT, read.PNEXT, read.TLEN,
                                  read.seq,
                                  read.raw_base_quality,
                                  phasing_info])
            save_file_fp.stdin.write(read_str + '\n')
            del aligned_reads[read_name]
        if read_fn != 'PIPE':
            save_file_fp.stdin.close()
            save_file_fp.wait()
    samtools_view_process.stdout.close()
    samtools_view_process.wait()

    if test_pos:
        save_file_fp = subprocess_popen(shlex.split("samtools index {}".format(
            read_fn + ('.{}_{}'.format(ctg_start, ctg_end) if is_ctg_range_given and not test_pos else ""))),
            stdin=PIPE, stdout=PIPE)
        save_file_fp.stdin.close()
        save_file_fp.wait()


def main():
    parser = ArgumentParser(description="Reads realignment")

    parser.add_argument('--bam_fn', type=str, default=None, required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--read_fn', type=str, default="PIPE",
                        help="Output realigned BAM. Default directly pass reads to CreateTensor_phasing using PIPE. Default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help="Realign reads only in the provided bed regions")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    # options for advanced users
    parser.add_argument('--minCoverage', type=float, default=2,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="EXPERIMENTAL: Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default: %(default)d")

    # options for debug purpose
    parser.add_argument('--extend_bed', type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    # options for internal process control
    ## Test in specific candidate position. Only for testing
    parser.add_argument('--test_pos', type=int, default=0,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    reads_realignment(args)


if __name__ == "__main__":
    main()
