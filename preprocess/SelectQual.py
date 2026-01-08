from sys import stdin
from argparse import ArgumentParser, SUPPRESS
import os

from shared.utils import file_path_from, log_warning

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

def select_phase_qual_from_stdin(args):

    """
    Select a global quality cut-off for phasing and reads haplotag.
    """
    qual_fn = args.qual_fn if args.qual_fn is not None else "phase_qual"
    var_pct_full = args.var_pct_full
    var_pct_phasing = args.var_pct_phasing
    low_qual_hete_var_pct = 1 - var_pct_phasing if var_pct_phasing is not None else var_pct_full
    phase_qual_list = []
    for row in stdin:
        if row[0] == '#':
            continue
        row = row.rstrip().split()
        ref_base, alt_base = row[3], row[4]
        # select heterozygous snp only
        if len(ref_base) != 1 or len(alt_base) != 1:
            continue
        qual, gt_info = row[5], row[9]
        genotype = gt_info.split(':')[0]
        if genotype == '0/1':
            phase_qual_list.append(float(qual))

    # in phase mode, var_pct_full is the proportion of low-quality heterozygous variants to be discarded for whatshap phasing
    phase_qual_list = sorted(phase_qual_list)
    low_phase_qual_list = phase_qual_list[:int(low_qual_hete_var_pct * len(phase_qual_list))]
    if len(low_phase_qual_list) == 0:
        print(log_warning(
            "[WARNING] Cannot find any 0/1 variant in pileup output using variant quality cut-off proportion: {}, total heterozygous variants: {}".format(
                low_qual_hete_var_pct, len(low_phase_qual_list))))
        print(log_warning("[WARNING] Set low variant quality score cut-off to 0.0"))
        qual_cut_off = 0.0
    else:
        qual_cut_off = low_phase_qual_list[-1]
    print ('[INFO] Select heterozygous pileup variants exceeding phasing quality cutoff {}'.format(round(qual_cut_off), 0))

    if args.output_fn:
        with open(os.path.join(args.output_fn, qual_fn), 'w') as output:
            output.write(str(qual_cut_off))



def select_qual_from_stdin(args):

    """
    Select a global quality cut-off for full alignment calling from pileup vcf file. False positive pileup variants
    and true variants missed by pileup calling would mostly have low quality score (reference quality score for missing
    variants), so only use a proportion of low quality variants for full alignment while maintain high quality pileup
    output, as full alignment calling is substantially slower than pileup calling.
    """
    var_pct_full = args.var_pct_full
    qual_fn = args.qual_fn if args.qual_fn is not None else "qual"
    vcf_fn = file_path_from(args.vcf_fn)
    ref_pct_full = args.ref_pct_full if args.ref_pct_full else var_pct_full
    # for efficiency, we use a maximum 30% reference candidates proportion for full-alignment calling, which is almost cover all false negative candidates
    # for ont platform, we set a default 10% reference candidates proportion for full-alignment calling unless a known vcf file is provided (genotyping mode)
    # directly set default value in run_clair3.sh from v0.1-r5
    # ref_pct_full = 0.1 if args.platform == 'ont' else ref_pct_full
    # ref_pct_full = min(ref_pct_full, 0.3)

    variant_qual_list = []
    ref_qual_list = []
    for row in stdin:
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
    low_variant_qual_list = variant_qual_list[:int(var_pct_full * len(variant_qual_list))]
    if len(low_variant_qual_list) == 0:
        print(log_warning(
            "[WARNING] Cannot find any low-quality 0/1 or 1/1 variant in pileup output using variant quality cut-off proportion: {}, total variants: {}".format(
                var_pct_full, len(variant_qual_list))))
        print(log_warning("[WARNING] Set low variant quality score cut-off to 0.0"))
        var_qual_cut_off = 0.0
    else:
        var_qual_cut_off = low_variant_qual_list[-1]

    # If a known vcf file is provided, use user-defined proportion
    low_ref_qual_list = ref_qual_list[:int(ref_pct_full * len(ref_qual_list))] if vcf_fn is None else ref_qual_list[:int(args.ref_pct_full * len(ref_qual_list))]
    if len(low_ref_qual_list) == 0:
        print(log_warning(
            "[WARNING] Cannot find any low-quality 0/0 reference calls in pileup output using reference quality cut-off proportion: {}, total reference calls: {}".format(
                ref_pct_full, len(ref_qual_list))))
        print(log_warning("[WARNING] Set low reference quality score cut-off to 0.0"))
        ref_qual_cut_off = 0.0
    else:
        ref_qual_cut_off = low_ref_qual_list[-1]
    print ('[INFO] Set variants quality cutoff {}'.format(round(var_qual_cut_off, 0)))
    print ('[INFO] Set reference calls quality cutoff {}'.format(round(ref_qual_cut_off, 0)))

    if args.output_fn:
        with open(os.path.join(args.output_fn, qual_fn), 'w') as output:
            output.write(str(var_qual_cut_off) + ' ' + str(ref_qual_cut_off))


def main():
    parser = ArgumentParser(description="Select quality cut-off for phasing and full alignment")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--output_fn', type=str, default=None, required=True,
                        help="Define the output folder, required")

    parser.add_argument('--var_pct_full', type=float, default=0.3,
                        help="Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3")

    parser.add_argument('--ref_pct_full', type=float, default=0.3,
                        help="Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont")

    parser.add_argument('--var_pct_phasing', type=float, default=0.7,
                        help="Specify an expected percentage of high quality 0/1 variants used in WhatsHap phasing, default: 0.8 for ont guppy5 and 0.7 for other platforms")

    parser.add_argument('--phase', action='store_true',
                        help="Select only heterozygous candidates for phasing or not, default: False")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    # options for internal process control
    ## Input the file that contains the quality cut-off for selecting low-quality pileup calls for phasing and full-alignment calling
    parser.add_argument('--qual_fn', type=str, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()
    if args.phase:
        select_phase_qual_from_stdin(args)
    else:
        select_qual_from_stdin(args)


if __name__ == "__main__":
    main()
