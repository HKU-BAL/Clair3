import shlex
import os
import sys
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict
from shared.intervaltree.intervaltree import IntervalTree


import shared.param_f as param
from shared.utils import subprocess_popen

def FiterHeteSnpPhasing(args):

    """
    Filter heterozygous snp variant for phasing, currently, we only filter snp variant with low quality socore as low
    quality variant contains more false positive variant that would lead to a larger minimum error correction loss.
    """
    qual_fn = args.qual_fn if args.qual_fn is not None else 'phase_qual'
    vcf_fn = args.vcf_fn
    var_pct_full = args.var_pct_full
    contig_name = args.ctgName
    split_folder = args.split_folder
    variant_dict = defaultdict(str)
    qual_set = defaultdict(int)
    found_qual_cut_off = False
    header = []

    #try to find the global quality cut off:
    f_qual = os.path.join(split_folder, qual_fn)
    if os.path.exists(f_qual):
        phase_qual_cut_off = float(open(f_qual, 'r').read().rstrip())
        found_qual_cut_off = True

    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
    for row in unzip_process.stdout:
        row = row.rstrip()
        if row[0] == '#':
            header.append(row + '\n')
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        ref_base = columns[3]
        alt_base = columns[4]
        genotype = columns[9].split(':')[0].replace('|', '/')

        if len(ref_base) == 1 and len(alt_base) == 1:
            if genotype == '0/1' or genotype=='1/0':
                variant_dict[pos] = row
                qual = float(columns[5])
                qual_set[pos] = qual

    if found_qual_cut_off:
        remove_low_qual_list = [[k,v] for k,v in qual_set.items() if v < phase_qual_cut_off ]
    else:
        remove_low_qual_list = sorted(qual_set.items(), key=lambda x: x[1])[:int(var_pct_full * len(qual_set))]
    for pos, qual in remove_low_qual_list:
        del variant_dict[pos]

    print ('[INFO] Total heterozygous SNP positions selected: {}: {}'.format(contig_name, len(variant_dict)))

    f = open(os.path.join(split_folder, '{}.vcf'.format(contig_name)), 'w')
    f.write(''.join(header))
    for key,row in sorted(variant_dict.items(), key=lambda x: x[0]):
        f.write(row +'\n')
    f.close()


def FiterHeteSnp_FP(args):

    """
    Filter heterozygous snp variant for calling, this is a testing function to validate various proportion of phasing
    effect on full alignment calling, currently for testing only.
    """

    vcf_fn = args.vcf_fn
    proportion = args.proportion
    chr_prefix = args.chr_prefix
    contig_name = args.ctgName
    phasing_window_size = param.phasing_window_size
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
    split_bed_size = args.split_bed_size
    split_folder = args.split_folder
    output = []
    snp = []
    need_phasing_list = []
    chr_prefix_length = len(chr_prefix)
    variant_dict = defaultdict(str)
    for row in unzip_process.stdout:

        if row[0] == '#':
            output.append(row.rstrip())
            continue
        columns = row.strip().split()

        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        ref_base = columns[3]
        alt_base = columns[4]
        genotype = columns[9].split(':')[0].replace('|', '/')
        qual = int(columns[5])
        if len(ref_base) == 1 and len(alt_base) == 1:
            if genotype == '0/1':
                snp.append((qual, pos))
                variant_dict[pos] = ref_base + '-' + alt_base
        else:
            need_phasing_list.append((qual, pos))

    qual_list = sorted(snp, key=lambda x: x[0])
    print('[INFO] Total hete snp variants:', len(qual_list))
    cut_off_index = int(len(qual_list) * proportion)
    hete_snp_row_list = [item[1] for item in qual_list[cut_off_index:]]
    print ('[INFO] Total hete snp filter matches:', len(hete_snp_row_list))

    qual_list = sorted(need_phasing_list, key=lambda x: -x[0])
    cut_off_index = int(len(qual_list) * proportion)
    need_phasing_row_list = sorted([item[1] for item in qual_list[cut_off_index:]])
    print('[INFO] Total variants need to be phased:', len(need_phasing_row_list))
    phasing_tree = IntervalTree()
    for item_idx, item in enumerate(need_phasing_list):
        pos = item[1]
        start = pos - phasing_window_size
        end = pos + phasing_window_size
        phasing_tree.addi(start, end)

    snp_tree = IntervalTree()
    for item in hete_snp_row_list:
        if len(phasing_tree.at(item)): snp_tree.addi(item, item + 1)

    region_num = len(need_phasing_row_list) // split_bed_size  + 1 if len(need_phasing_row_list) % split_bed_size else len(need_phasing_row_list) // split_bed_size
    for idx in range(region_num):
        split_output = need_phasing_row_list[idx * split_bed_size : (idx+1) * split_bed_size]
        start, end = split_output[0] - phasing_window_size, split_output[-1] + phasing_window_size
        overlaps = snp_tree.overlap(start, end)
        snp_split_out = []
        for overlap in overlaps:
            snp_split_out.append((overlap[0], overlap[0] + 1, 1))
        split_output = [(item - param.flankingBaseNum, item+1 + param.flankingBaseNum, 0) for item in split_output] # a windows region for create tensor
        print (len(split_output), len(snp_split_out))
        split_output += snp_split_out
        split_output = sorted(split_output, key=lambda x: x[0])

        with open(os.path.join(split_folder, 'split_{}.{}'.format(contig_name[chr_prefix_length:], idx)), 'w') as output_file:
            output_file.write('\n'.join(['\t'.join([contig_name, str(x[0]-1), str(x[1]-1), str(x[2]), variant_dict[x[0]]]) for x in split_output]) + '\n') # bed format

def FiterHeteSnp(args):

    """
    Filter heterozygous snp variant for training, if there are too many candidates for full alignment training, we
    would select more in low quality variants, which is more challenging for pileup model to predict and using more
    information will benefit calling those variants.
    """

    vcf_fn = args.vcf_fn # true vcf var
    alt_fn = args.alt_fn
    var_pct_full = args.var_pct_full
    ref_pct_full = args.ref_pct_full if args.ref_pct_full is not None else var_pct_full
    chr_prefix = args.chr_prefix
    contig_name = args.ctgName
    phasing_window_size = param.phasing_window_size
    chunk_id = args.chunk_id - 1 if args.chunk_id else None # 1-base to 0-base
    DEPTH = args.depth
    chunk_num = args.chunk_num
    sample_name = args.sampleName
    split_bed_size = args.split_bed_size
    split_folder = args.split_folder
    extend_bp = param.extend_bp
    phasing_info_in_bam = args.phasing_info_in_bam
    need_phasing_list = []
    need_phasing_set = set()
    ref_call_pos_list = []
    chr_prefix_length = len(chr_prefix)
    variant_dict = defaultdict(str)
    realign_window_size = args.realign_window_size if args.realign_window_size is not None else param.flankingBaseNum
    candidate_positions = set()

    if vcf_fn and os.path.exists(vcf_fn):
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
        for row in unzip_process.stdout:
            row =row.rstrip()
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
            if '|' not in genotype: #unphasable
                continue
            variant_dict[pos] = '-'.join([ref_base, alt_base, ('1' if genotype == '0|1' else '2'), phase_set])

    if alt_fn is not None:
        # vcf format
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (alt_fn)))
        for row in unzip_process.stdout:
            if row[0] == '#':
                continue
            columns =row.rstrip().split('\t')
            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            qual = float(columns[5])


            candidate_positions.add(pos)
            #ref_call was marked as '.' after v0.1-r5
            if ref_base == alt_base or alt_base == ".":
                ref_call_pos_list.append((pos,qual))
            else:
                need_phasing_list.append((pos,qual))
                need_phasing_set.add(pos)

        low_qual_ref_list = sorted(ref_call_pos_list, key=lambda x: x[1])[:int(ref_pct_full * len(ref_call_pos_list))]
        low_qual_variant_list = sorted(need_phasing_list, key=lambda x: x[1])[:int(var_pct_full * len(need_phasing_list))]

        #calling with phasing_info_in_bam: select low qual ref and low qual vairant for phasing calling
        if phasing_info_in_bam:
            print('[INFO] {} {} total low qual ref calling to process: {}'.format(sample_name, contig_name, len(low_qual_ref_list)))
            print('[INFO] {} {} total low qual variant calling to process: {}'.format(sample_name, contig_name, len(low_qual_variant_list)))

            need_phasing_row_list = set([item[0] for item in low_qual_ref_list] + [item[0] for item in low_qual_variant_list])
            need_phasing_row_list = sorted(list(need_phasing_row_list))

            if chunk_num:
                all_candidate_size = len(need_phasing_row_list)
                chunk_size = all_candidate_size // chunk_num + 1 if all_candidate_size % chunk_num else all_candidate_size // chunk_num

                for chunk_idx in range(chunk_num):
                    start_pos = chunk_idx * chunk_size
                    end_pos = min(start_pos + chunk_size, all_candidate_size)
                    split_output = need_phasing_row_list[start_pos:end_pos]
                    split_output = [(item - realign_window_size, item + realign_window_size + 2) for item in
                                    split_output]  # a windows region for create tensor # samtools mpileup not include last position

                    split_output = sorted(split_output, key=lambda x: x[0])
                    with open(os.path.join(split_folder,
                                           '{}_{}_{}_{}'.format(sample_name, DEPTH, contig_name[chr_prefix_length:], chunk_idx+1)), # zero-base to one-base
                              'w') as output_file:
                        output_file.write('\n'.join(
                            ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1), ]) for x in
                             split_output]) + '\n')  # bed format
                return

            region_num = len(need_phasing_row_list) // split_bed_size + 1 if len(
                need_phasing_row_list) % split_bed_size else len(need_phasing_row_list) // split_bed_size
            for idx in range(region_num):
                split_output = need_phasing_row_list[idx * split_bed_size: (idx + 1) * split_bed_size]
                split_output = [(item - realign_window_size, item + realign_window_size + 2) for item in
                                split_output]  # a windows region for create tensor # samtools mpileup not include last position

                split_output = sorted(split_output, key=lambda x: x[0])

                with open(os.path.join(split_folder, '{}.{}_{}'.format(contig_name[chr_prefix_length:], split_output[0][0], split_output[-1][1])),
                          'w') as output_file:
                    output_file.write('\n'.join(
                        ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1),]) for x in
                         split_output]) + '\n')  # bed format
            return

        for pos, qual in low_qual_ref_list:
            need_phasing_set.add(pos)

    # train or call in all_pos
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

    region_num = len(need_phasing_row_list) // split_bed_size + 1 if len(need_phasing_row_list) % split_bed_size else len(need_phasing_row_list) // split_bed_size
    for idx in range(region_num):
        split_output = need_phasing_row_list[idx * split_bed_size : (idx+1) * split_bed_size]

        start = split_output[0]
        end = split_output[-1]
        extend_start, extend_end = start - phasing_window_size, end + phasing_window_size
        overlaps = snp_tree.overlap(extend_start, extend_end)
        snp_split_out = []
        for overlap in overlaps:
            snp_split_out.append((contig_name, overlap[0] - extend_bp - 1 - 1, overlap[0] + 1 + extend_bp - 1, variant_dict[overlap[0]]))# bed format
        split_output = [(contig_name, item - realign_window_size-1, item+realign_window_size+1-1) for item in split_output] # a windows region for create tensor # bed format

        split_output += snp_split_out
        split_output = sorted(split_output, key=lambda x: x[1])

        with open(os.path.join(split_folder, '{}.{}_{}'.format(contig_name[chr_prefix_length:], start, end)), 'w') as output_file:
            output_file.write('\n'.join(['\t'.join(map(str, x)) for x in split_output]) + '\n') # bed format

def main():
    parser = ArgumentParser(description="Select heterozygous snp candidates for WhatsHap phasing")

    parser.add_argument('--split_folder', type=str, default=None,
                        help="Path to directory that stores small bed region for raw alignment. (default: %(default)s)")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Path of the input vcf file. (default: %(default)s)")

    parser.add_argument('--var_pct_full', type=float, default=0.3,
                        help="Default variant call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ref_pct_full', type=float, default=None,
                        help="Default reference call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--phase', action='store_false',
                        help="Only select hete candidates for phasing, default: True")

    parser.add_argument('--sampleName', type=str, default="",
                        help="Define the sample name to be shown in the VCF file, optional")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Input bam or sam have phasing info in HP tag, default: False")

    parser.add_argument('--split_bed_size', type=int, default=1000,
                        help="DEBUG: Default split bed size for parallel excution, default: %(default)s")

    parser.add_argument('--calling', type=int, default=0,
                        help="DEBUG: Path of the output folder, default: %(default)s")

    parser.add_argument('--realign_window_size', type=int, default=None,
                        help="DEBUG: The window size of read realignment, work with need_realignment")

    parser.add_argument('--split_region_size', type=int, default=40000000,
                        help="DEBUG: Vcf phasing split_region_size default: %(default)s")

    # options for internal process control
    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Output all alternative candidates path
    parser.add_argument('--all_alt_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Default chr prefix for contig name
    parser.add_argument('--chr_prefix', type=str, default='chr',
                        help=SUPPRESS)

    ## Default subsample depth for subsample bam file, 1000 means no subsampling
    parser.add_argument('--depth', type=int, default=1000,
                        help=SUPPRESS)

    ## Path of provided alternative file
    parser.add_argument('--alt_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Input the file that contains the quality cut-off for selecting low-quality pileup calls for phasing and full-alignment calling
    parser.add_argument('--qual_fn', type=str, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    #
    if args.phase:
        FiterHeteSnpPhasing(args)
    elif args.calling == 1:
        FiterHeteSnp_FP(args)
    else:
        FiterHeteSnp(args)

if __name__ == "__main__":
    main()
