from sys import stdin
from argparse import ArgumentParser
import os

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def select_phase_qual_from_stdin(args):

    """
    Select a global quality cut-off for phasing and reads haplotag.
    """

    phase_qual_list = []
    for row in stdin.readlines():
        if row[0] == '#':
            continue
        row = row.rstrip().split()
        ref_base, alt_base = row[3], row[4]
        # select hete snp only
        if len(ref_base) != 1 or len(alt_base) != 1:
            continue
        qual, gt_info = row[5], row[9]
        genotype = gt_info.split(':')[0]
        if genotype == '0/1':
            phase_qual_list.append(float(qual))

    phase_qual_list = sorted(phase_qual_list)
    qual_cut_off = phase_qual_list[:int(args.var_pct_full * len(phase_qual_list))][-1]
    print ('[INFO] Select phasing quality cut off {}'.format(round(qual_cut_off), 0))

    if args.output_fn:
        with open(os.path.join(args.output_fn, 'phase_qual'), 'w') as output:
            output.write(str(qual_cut_off))



def select_qual_from_stdin(args):

    """
    Select a global quality cut-off for full alignment directly from pileup vcf file. False positive pileup variants
    and true variants missed by pileup calling would mostly have low quality score (reference quality score for missing
    variants), so only use a proportion of low quality variants for full alignment while maintain high quality pileup
    output, as full alignment calling is substantially slower than pileup calling.
    """
    var_pct_full = args.var_pct_full
    ref_pct_full = args.ref_pct_full if args.ref_pct_full else var_pct_full
    ref_pct_full = 0.1 if  args.platform == 'ont' else ref_pct_full
    variant_qual_list = []
    ref_qual_list = []
    for row in stdin.readlines():
        if row[0] == '#':
            continue
        row = row.rstrip().split()

        qual, gt_info = row[5], row[9]
        genotype = gt_info.split(':')[0]
        if genotype == '0/0':
            ref_qual_list.append(float(qual))
        else:
            variant_qual_list.append(float(qual))

    ref_qual_list = sorted(ref_qual_list)
    variant_qual_list = sorted(variant_qual_list)
    var_qual_cut_off = variant_qual_list[:int(var_pct_full * len(variant_qual_list))][-1]
    #for efficency, we use maxinum 30% reference qual, which is almost cover all false negative candidates
    ref_qual_cut_off = ref_qual_list[:int(min(ref_pct_full, 0.3) * len(ref_qual_list))][-1]
    print ('[INFO] Select variant quality cut off {}'.format(round(var_qual_cut_off, 0)))
    print ('[INFO] Select reference quality cut off {}'.format(round(ref_qual_cut_off, 0)))

    if args.output_fn:
        with open(os.path.join(args.output_fn, 'qual'), 'w') as output:
            output.write(str(var_qual_cut_off) + ' ' + str(ref_qual_cut_off))


def main():
    parser = ArgumentParser(description="Select quality cut-off for phasing and full alignment")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--output_fn', type=str, default=None, required=True,
                        help="Define the output folder, required")

    parser.add_argument('--var_pct_full', type=float, default=0.3,
                        help="Default variant call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ref_pct_full', type=float, default=None,
                        help="Default reference call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--phase', action='store_true',
                        help="Select only heterozygous candidates for phasing or not, default: False")

    args = parser.parse_args()
    if args.phase:
        select_phase_qual_from_stdin(args)
    else:
        select_qual_from_stdin(args)


if __name__ == "__main__":
    main()
