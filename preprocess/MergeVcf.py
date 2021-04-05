import sys
import os
import shlex
from argparse import ArgumentParser
import logging

logging.basicConfig(format='%(message)s', level=logging.INFO)

from shared.utils import subprocess_popen, str2bool
import shared.param_f as param
from shared.interval_tree import bed_tree_from, is_region_in
from clair3.MakeGvcf import gvcfGenerator


def filtration_value_from(quality_score_for_pass, quality_score):
    if quality_score_for_pass is None:
        return "."
    if quality_score <= quality_score_for_pass:
        return "LowQual"
    return "PASS"


def MarkLowQual(row, quality_score_for_pass, qual):
    row = row.split("\t")
    row[6] = filtration_value_from(quality_score_for_pass, qual)
    return '\t'.join(row)


def MergeVcf_illumina(args):
    # region vcf merge for illumina, as read realignment will make candidate varaints shift and missing.
    candidate_bed_fn = args.candidate_bed_fn
    output_fn = args.output_fn
    full_alignment_vcf_fn = args.full_alignment_vcf_fn
    pileup_vcf_fn = args.pileup_vcf_fn  # true vcf var
    contig_name = args.ctgName

    tree = bed_tree_from(bed_file_path=candidate_bed_fn, padding=param.no_of_positions)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (pileup_vcf_fn)))
    output = []
    pileup_count = 0
    for row in unzip_process.stdout:
        if row[0] == '#':
            output.append(row)
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name != None and ctg_name != contig_name:
            continue
        pos = int(columns[1])
        pass_bed = is_region_in(tree, ctg_name, pos)
        if not pass_bed:
            output.append(row)
            pileup_count += 1

    unzip_process.stdout.close()
    unzip_process.wait()

    realigned_vcf_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (full_alignment_vcf_fn)))
    realiged_read_num = 0
    for row in realigned_vcf_unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name != None and ctg_name != contig_name:
            continue

        pos = int(columns[1])
        if is_region_in(tree, ctg_name, pos):
            output.append(row)
            realiged_read_num += 1

    logging.info('[INFO] Pileup positions varaints proceeded in {}: {}'.format(contig_name, pileup_count))
    logging.info('[INFO] Realiged positions varaints proceeded in {}: {}'.format(contig_name, realiged_read_num))
    realigned_vcf_unzip_process.stdout.close()
    realigned_vcf_unzip_process.wait()

    with open(output_fn, 'w') as output_file:
        output_file.write(''.join(output))


def MergeVcf(args):
    """
    Merge pileup and full alignment vcf output. As current workflow only use part of low quality score candidate
    variants of pileup output for full alignment.
    """

    output_fn = args.output_fn
    full_alignment_vcf_fn = args.full_alignment_vcf_fn
    pileup_vcf_fn = args.pileup_vcf_fn  # true vcf var
    contig_name = args.ctgName
    QUAL = args.qual
    full_alignment_vcf_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (full_alignment_vcf_fn)))

    pileup_output = []
    full_alignment_output = []
    full_alignment_output_set = set()

    for row in full_alignment_vcf_unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name != None and ctg_name != contig_name:
            continue
        pos = int(columns[1])
        qual = int(float(columns[5]))
        row = MarkLowQual(row, QUAL, qual)
        full_alignment_output.append(((ctg_name, pos), row, qual))

        full_alignment_output_set.add((ctg_name, pos))

    full_alignment_vcf_unzip_process.stdout.close()
    full_alignment_vcf_unzip_process.wait()

    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (pileup_vcf_fn)))

    header = []

    for row in unzip_process.stdout:
        if row[0] == '#':
            header.append(row)
            continue

        columns = row.rstrip().split('\t')
        ctg_name = columns[0]
        if contig_name and contig_name != ctg_name:
            continue
        pos = int(columns[1])
        qual = int(float(columns[5]))

        if (ctg_name, pos) in full_alignment_output_set:
            continue
        row = MarkLowQual(row, QUAL, qual)
        pileup_output.append(((ctg_name, pos), row))

    logging.info('[INFO] Pileup varaints proceeded in {}: {}'.format(contig_name, len(pileup_output)))
    logging.info('[INFO] Full alignemnt varaints proceeded in {}: {}'.format(contig_name, len(full_alignment_output)))

    with open(output_fn, 'w') as output_file:
        output_file.write(''.join(header))
        output = [item[1] for item in sorted(pileup_output + full_alignment_output, key=lambda x: x[0])]
        output_file.write(''.join(output))


def mergeNonVariant(args):
    '''
    merge the variant calls and non-variants

    '''
    gvcf_generator = gvcfGenerator(ref_path=args.ref_fn, samtools=args.samtools)
    raw_gvcf_path = args.non_var_gvcf_fn
    raw_vcf_path = args.output_fn
    if (args.gvcf_fn == None):
        save_path = args.call_fn.split('.')[0] + '.g.vcf'
    else:
        save_path = args.gvcf_fn
        logging.info("merge calls")
        gvcf_generator.mergeCalls(raw_vcf_path, raw_gvcf_path, save_path, args.sampleName, args.ctgName, args.ctgStart,
                                  args.ctgEnd)

    # remove temporary file
    # cmd = "rm " + raw_gvcf_path
    # os.system(cmd)

    pass


def main():
    logging.info("Merge pileup and full alignment vcf")
    parser = ArgumentParser(description="Generate 1-based variant candidates using alignments")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Path to directory that stores small bins, default: %(default)s")

    parser.add_argument('--candidate_bed_fn', type=str, default=None,
                        help="Path of the output folder. default: %(default)s")

    parser.add_argument('--full_alignment_vcf_fn', type=str, default=None,
                        help="Path of full alignment vcf file, default: %(default)s")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Path of pileup vcf file, default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--qual', type=int, default=0,
                        help="If set, variants with ≥$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--platform', type=str, default=None,
                        help="Sequencing platform of the input. Options: 'ont,pb,illumina', default: %(default)s")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    parser.add_argument('--non_var_gvcf_fn', type=str, default=None,
                        help='path of non-variant gvcf')

    parser.add_argument('--gvcf_fn', type=str, default=None,
                        help='File name of gvcf outpupt')

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, optional")

    parser.add_argument('--samtools', type=str, default='samtools',
                        help="Path to the 'samtools', samtools verision >= 1.10 is required, default: %(default)s")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file")

    # options for debug purpose
    parser.add_argument('--ctgStart', type=int, default=None,
                        help="DEBUG: The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="DEBUG: The 1-based inclusive ending position of the sequence to be processed")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    # realignemnt region merge
    if args.platform == 'illumina':
        MergeVcf_illumina(args=args)
    else:
        MergeVcf(args=args)
    if (args.gvcf):
        mergeNonVariant(args)


if __name__ == "__main__":
    main()
