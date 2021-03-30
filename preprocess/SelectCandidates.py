import shlex
from argparse import ArgumentParser
import shared.param_f as param
from collections import defaultdict
import os
from intervaltree import IntervalTree
import math
import sys
import logging

logging.basicConfig(format='%(message)s', level=logging.INFO)

from shared.utils import subprocess_popen, IUPAC_base_to_num_dict as BASE2NUM, region_from, reference_sequence_from, str2bool

def gaussian_distribution(x, mu, sig=16):
    return math.exp(-math.pow(x - mu, 2.) / (2 * math.pow(sig, 2.)))


def discrete_gaussian_pro(entropy_windnow):
    gaussian_pro = [gaussian_distribution(index, entropy_windnow / 2, 1) for index in range(entropy_windnow)]
    return gaussian_pro


def calculate_sequence_entropy(sequence, entropy_window=None, kmer=5):
    """
    We use a kmer-based sequence entropy calculation to measure the complexity of a region.
    sequence: a chunked sequence around a candidate position, default no_of_positions = flankingBaseNum + 1 + flankingBaseNum
    entropy_window: a maximum entropy window for scanning, if the sequence is larger than the entropy window, a slide
    window would be adopted for measurement.
    kmer: default kmer size for sequence entropy calculation.
    """

    count_repeat_kmer_counts = [0] * (entropy_window + 2)
    count_repeat_kmer_counts[0] = entropy_window

    entropy = [0.0] * (entropy_window + 2)
    for i in range(1, entropy_window + 2):
        e = 1.0 / entropy_window * i
        entropy[i] = e * math.log(e)
    entropy_mul = -1 / math.log(entropy_window)
    entropy_kmer_space = 1 << (2 * kmer)

    kmer_hash_counts = [0] * entropy_kmer_space  # value should smaller than len(seq)
    mask = -1 if kmer > 15 else ~((-1) << (2 * kmer))
    kmer_suffix, kmer_prefix = 0, 0

    i = 0
    i2 = -entropy_window
    entropy_sum = 0.0
    all_entropy_sum = [0.0] * len(sequence)
    while (i2 < len(sequence)):

        if (i < len(sequence)):
            n = BASE2NUM[sequence[i]]
            kmer_suffix = ((kmer_suffix << 2) | n) & mask

            count_repeat_kmer_counts[kmer_hash_counts[kmer_suffix]] -= 1
            entropy_sum -= entropy[kmer_hash_counts[kmer_suffix]]
            kmer_hash_counts[kmer_suffix] += 1
            count_repeat_kmer_counts[kmer_hash_counts[kmer_suffix]] += 1
            entropy_sum += entropy[kmer_hash_counts[kmer_suffix]]

        if i2 >= 0 and i < len(sequence):
            n2 = BASE2NUM[sequence[i2]]
            kmer_prefix = ((kmer_prefix << 2) | n2) & mask  # add base info
            count_repeat_kmer_counts[kmer_hash_counts[kmer_prefix]] -= 1
            entropy_sum -= entropy[kmer_hash_counts[kmer_prefix]]
            kmer_hash_counts[kmer_prefix] -= 1
            count_repeat_kmer_counts[kmer_hash_counts[kmer_prefix]] += 1
            entropy_sum += entropy[kmer_hash_counts[kmer_prefix]]
            all_entropy_sum[i] = entropy_sum
        i += 1
        i2 += 1
    return entropy_sum * entropy_mul


def sqeuence_entropy_from(samtools_execute_command, fasta_file_path, contig_name, candidate_positions):
    """
    Calculate sequence entropy in a specific candidate windows, variants in low sequence entropy regions (low
    mappability regions, such as homopolymer, tandem repeat, segmental duplications regions) would more likely have
    more complex variants representation, which is beyond pileup calling. Hence, those candidate variants are re-called by
    full alignment calling.
    We use a kmer-based sequence entropy calculation to measure the complexity of a region, we would directly query the
    chunked reference sequence for sequence entropy calculation for each candidate variant.
    """

    ref_regions = []
    reference_start, reference_end = min(list(candidate_positions)) - param.no_of_positions, max(
        list(candidate_positions)) + param.no_of_positions + 1
    reference_start = 1 if reference_start < 1 else reference_start
    ref_regions.append(region_from(ctg_name=contig_name, ctg_start=reference_start, ctg_end=reference_end))
    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools_execute_command,
        fasta_file_path=fasta_file_path,
        regions=ref_regions
    )
    if reference_sequence is None or len(reference_sequence) == 0:
        sys.exit("[ERROR] Failed to load reference seqeunce from file ({}).".format(fasta_file_path))

    entropy_window = param.no_of_positions
    candidate_positions_entropy_list = []
    for pos in candidate_positions:
        ref_seq = reference_sequence[
                  pos - param.flankingBaseNum - reference_start: pos + param.flankingBaseNum + 1 - reference_start]
        sequence_entropy = calculate_sequence_entropy(sequence=ref_seq, entropy_window=entropy_window)
        candidate_positions_entropy_list.append((pos, sequence_entropy))

    return candidate_positions_entropy_list


def SelectCandidates(args):
    """
    Select low quality and low sequence entropy candidate variants for full aligement. False positive pileup variants
    and true variants missed by pileup calling would mostly have low quality score (reference quality score for missing
    variants), so only use a proportion of low quality variants for full alignment while maintain high quality pileup
    output, as full alignment calling is substantially slower than pileup calling.
    """

    vcf_fn = args.vcf_fn  # true vcf var
    alt_fn = args.alt_fn
    proportion = args.pro
    seq_entropy_pro = args.seq_entropy_pro
    contig_name = args.ctgName
    phasing_window_size = param.phasing_window_size
    platform = args.platform
    split_bed_size = args.split_bed_size
    split_folder = args.split_folder
    extend_bp = param.extend_bp
    call_low_seq_entropy = args.call_low_seq_entropy
    phasing_info_in_bam = args.phasing_info_in_bam
    need_phasing_list = []
    need_phasing_set = set()
    ref_call_pos_list = []
    variant_dict = defaultdict(str)
    realign_window_size = args.realign_window_size if args.realign_window_size is not None else param.flankingBaseNum

    fasta_file_path = args.ref_fn
    samtools_execute_command = args.samtools

    found_qual_cut_off = False
    low_sequence_entropy_list = []
    # try to find the global quality cut off:
    f_qual = os.path.join(split_folder, 'qual')
    if os.path.exists(f_qual):
        line = open(f_qual, 'r').read().rstrip().split(' ')
        var_qual, ref_qual = float(line[0]), float(line[1])
        found_qual_cut_off = True

    if vcf_fn and os.path.exists(vcf_fn):
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
        for row in unzip_process.stdout:
            row = row.rstrip()
            if row[0] == '#':
                continue
            columns = row.strip().split('\t')

            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            genotype_info = columns[9].split(':')
            genotype, phase_set = genotype_info[0], genotype_info[-1]
            if '|' not in genotype:  # unphasable
                continue
            variant_dict[pos] = '-'.join([ref_base, alt_base, ('1' if genotype == '0|1' else '2'), phase_set])

    if alt_fn is not None:
        # vcf format
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (alt_fn)))
        for row in unzip_process.stdout:
            if row[0] == '#':
                continue
            columns = row.rstrip().split('\t')
            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            qual = float(columns[5])

            # reference calling
            if ref_base == alt_base:
                ref_call_pos_list.append((pos, qual))
            else:
                need_phasing_list.append((pos, qual))
                need_phasing_set.add(pos)

        if found_qual_cut_off:
            low_qual_ref_list = [[k, v] for k, v in ref_call_pos_list if v < ref_qual]
            low_qual_variant_list = [[k, v] for k, v in need_phasing_list if v < var_qual]
        else:
            low_qual_ref_list = sorted(ref_call_pos_list, key=lambda x: x[1])[:int(proportion * len(ref_call_pos_list))]
            low_qual_variant_list = sorted(need_phasing_list, key=lambda x: x[1])[
                                    :int(proportion * len(need_phasing_list))]
        
        if call_low_seq_entropy:
            candidate_positions = sorted(ref_call_pos_list, key=lambda x: x[1])[
                                  :int((proportion + seq_entropy_pro) * len(ref_call_pos_list))] + sorted(need_phasing_list,
                                                                                                          key=lambda x: x[
                                                                                                              1])[:int(
                (proportion + seq_entropy_pro) * len(need_phasing_list))]
            candidate_positions = set([item[0] for item in candidate_positions])
    
            candidate_positions_entropy_list = sqeuence_entropy_from(samtools_execute_command=samtools_execute_command,
                                                                     fasta_file_path=fasta_file_path,
                                                                     contig_name=contig_name,
                                                                     candidate_positions=candidate_positions)
    
            low_sequence_entropy_list = sorted(candidate_positions_entropy_list, key=lambda x: x[1])[
                                        :int(seq_entropy_pro * len(candidate_positions_entropy_list))]

        # calling with phasing_info_in_bam: select low qual ref and low qual vairant for phasing calling
        if phasing_info_in_bam:
            logging.info(
                '[INFO] Total low qual ref variants to procceed in {}: {}'.format(contig_name, len(low_qual_ref_list)))
            logging.info(
                '[INFO] Total low qual variants to procceed in {}: {}'.format(contig_name, len(low_qual_variant_list)))
            if call_low_seq_entropy:
                logging.info('[INFO] Total low sequence entropy variants to procceed in {}: {}'.format(contig_name, len(
                    low_sequence_entropy_list)))

            need_phasing_row_list = set(
                [item[0] for item in low_qual_ref_list] + [item[0] for item in low_qual_variant_list] + [item[0] for
                                                                                                         item in
                                                                                                         low_sequence_entropy_list])
            need_phasing_row_list = sorted(list(need_phasing_row_list))

            region_num = len(need_phasing_row_list) // split_bed_size + 1 if len(
                need_phasing_row_list) % split_bed_size else len(need_phasing_row_list) // split_bed_size

            for idx in range(region_num):
                # a windows region for create tensor # samtools mpileup not include last position
                split_output = need_phasing_row_list[idx * split_bed_size: (idx + 1) * split_bed_size]

                if platform == 'illumina':
                    region_size = param.region_size
                    split_output = [(item // region_size * region_size - param.no_of_positions,
                                     item // region_size * region_size + region_size + param.no_of_positions) for item
                                    in split_output]
                else:
                    split_output = [(item - realign_window_size, item + realign_window_size + 2) for item in
                                    split_output]

                split_output = sorted(split_output, key=lambda x: x[0])

                # currently deprecate using ctgName.start_end as file name, which will run similar regions for several times when start and end has slight difference
                # output_path = os.path.join(split_folder, '{}.{}_{}'.format(contig_name, split_output[0][0], split_output[-1][1]))
                output_path = os.path.join(split_folder, '{}.{}_{}'.format(contig_name, idx, region_num))
                with open(output_path, 'w') as output_file:
                    output_file.write('\n'.join(
                        ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1), ]) for x in
                         split_output]) + '\n')  # bed format
            return

        for pos, qual in low_qual_ref_list:
            need_phasing_set.add(pos)

    # Call variant in all candidate position
    elif args.all_alt_fn is not None:
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (args.all_alt_fn)))
        for row in unzip_process.stdout:
            if row[0] == '#':
                continue
            columns = row.rstrip().split('\t')
            ctg_name, pos = columns[0].split()
            pos = int(pos)
            if contig_name and contig_name != ctg_name:
                continue
            need_phasing_set.add(pos)

    need_phasing_row_list = sorted(list(set(need_phasing_set)))
    snp_tree = IntervalTree()
    hete_snp_row_list = sorted(list(set(variant_dict.keys()).intersection(set(need_phasing_row_list))))
    print('[INFO] Total hete snp with reads support in {}: '.format(contig_name), len(hete_snp_row_list))
    print('[INFO] Total candidates need to be processed in {}: '.format(contig_name), len(need_phasing_row_list))

    for item in hete_snp_row_list:
        snp_tree.addi(item, item + 1)

    region_num = len(need_phasing_row_list) // split_bed_size + 1 if len(
        need_phasing_row_list) % split_bed_size else len(need_phasing_row_list) // split_bed_size
    for idx in range(region_num):
        split_output = need_phasing_row_list[idx * split_bed_size: (idx + 1) * split_bed_size]

        start = split_output[0]
        end = split_output[-1]
        extend_start, extend_end = start - phasing_window_size, end + phasing_window_size
        overlaps = snp_tree.overlap(extend_start, extend_end)
        snp_split_out = []
        for overlap in overlaps:
            snp_split_out.append((contig_name, overlap[0] - extend_bp - 1 - 1, overlap[0] + 1 + extend_bp - 1,
                                  variant_dict[overlap[0]]))  # bed format
        split_output = [(contig_name, item - realign_window_size - 1, item + realign_window_size + 1 - 1) for item in
                        split_output]  # a windows region for create tensor # bed format

        split_output += snp_split_out
        split_output = sorted(split_output, key=lambda x: x[1])

        with open(os.path.join(split_folder, '{}.{}_{}'.format(contig_name, start, end)), 'w') as output_file:
            output_file.write('\n'.join(['\t'.join(map(str, x)) for x in split_output]) + '\n')  # bed format


def main():
    parser = ArgumentParser(description="Select pileup candidates for full alignment")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Path to directory that stores small bins. (default: %(default)s)")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select specific platform for variant calling. Optional: 'ont,pb,illumina', default: %(default)s")

    parser.add_argument('--split_folder', type=str, default=None,
                        help="Path to directory that stores candidate region. (default: %(default)s)")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Path of the output folder. (default: %(default)s)")

    parser.add_argument('--pro', type=float, default=0.2,
                        help="Define the proportion for full alignment calling. (default: %(default)s)")
    
    parser.add_argument('--ref_pro', type=float, default=0.2,
                        help="Define the reference proportion for full alignment calling. (default: %(default)s)")

    parser.add_argument('--call_low_seq_entropy', type=str2bool, default=False,
                        help="Whether add indel length for training and calling, default true for raw alignment")

    parser.add_argument('--seq_entropy_pro', type=float, default=0.05,
                        help="Define the proportion for select low sequence entropy proportion for full alignment calling. (default: %(default)s)")

    parser.add_argument('--chr_prefix', type=str, default='chr',
                        help="Define ctgName prefix for processing. (default: %(default)s)")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--region_size', type=int, default=param.region_size,
                        help="Define the region size for illumina realignment calling. (default: %(default)s)")

    parser.add_argument('--split_bed_size', type=int, default=10000,
                        help="Define the candidate bed size for each split bed file. (default: %(default)s)")

    parser.add_argument('--alt_fn', type=str, default=None,
                        help="Path of pileup alternative vcf file. (default: %(default)s)")

    parser.add_argument('--realign_window_size', type=int, default=None,
                        help="Define the realign window size for long read phasing and realignment, only for testing. (default: %(default)s)")

    parser.add_argument('--phasing_info_in_bam', action='store_false',
                        help="Whether input bam or sam have phasing info in HP tag, default: False")

    parser.add_argument('--all_alt_fn', type=str, default=None,
                        help="Path of the all alternative vcf folder, only for testing. (default: %(default)s)")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools verision >= 1.10 is required, default: %(default)s")

    args = parser.parse_args()

    SelectCandidates(args)


if __name__ == "__main__":
    main()
