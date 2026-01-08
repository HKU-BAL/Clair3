import sys
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
from shared.utils import subprocess_popen, vcf_candidates_from

class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()

def OutputVariant(args):
    var_fn = args.var_fn
    vcf_fn = args.vcf_fn
    truth_vcf_fn = args.truth_vcf_fn
    ctg_name = args.ctgName
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd

    truth_vcf_set = set()
    variant_set = set()
    if args.truth_vcf_fn is not None:
        truth_vcf_set = set(vcf_candidates_from(vcf_fn=truth_vcf_fn, contig_name=ctg_name))
    if args.var_fn != "PIPE":
        var_fpo = open(var_fn, "wb")
        var_fp = subprocess_popen(shlex.split("gzip -c"), stdin=PIPE, stdout=var_fpo)
    else:
        var_fp = TruthStdout(sys.stdout)

    is_ctg_region_provided = ctg_start is not None and ctg_end is not None

    vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))

    for row in vcf_fp.stdout:
        columns = row.strip().split()
        if columns[0][0] == "#":
            continue

        # position in vcf is 1-based
        chromosome, position = columns[0], columns[1]
        if chromosome != ctg_name:
            continue
        if is_ctg_region_provided and not (ctg_start <= int(position) <= ctg_end):
            continue
        reference, alternate, last_column = columns[3], columns[4], columns[-1]
        # normal GetTruth
        genotype = last_column.split(":")[0].replace("/", "|").replace(".", "0").split("|")
        genotype_1, genotype_2 = genotype

        # 1000 Genome GetTruth (format problem) (no genotype is given)
        if int(genotype_1) > int(genotype_2):
            genotype_1, genotype_2 = genotype_2, genotype_1

        #remove * to guarentee vcf match
        if '*' in alternate:
            alternate = alternate.split(',')
            if int(genotype_1) + int(genotype_2) != 3 or len(alternate) != 2:
                print ('error with variant represatation')
                continue
            alternate = ''.join([alt_base for alt_base in alternate if alt_base != '*'])
            # * always have a genotype 1/2

            genotype_1, genotype_2 = '0', '1'

        variant_set.add(int(position))
        var_fp.stdin.write(" ".join((chromosome, position, reference, alternate, genotype_1, genotype_2)))
        var_fp.stdin.write("\n")

    for position in truth_vcf_set:
        if position not in variant_set:
            # miss variant set used in Tensor2Bin
            var_fp.stdin.write(" ".join((chromosome, str(position), "None", "None", "-1", "-1")))
            var_fp.stdin.write("\n")

    vcf_fp.stdout.close()
    vcf_fp.wait()

    if args.var_fn != "PIPE":
        var_fp.stdin.close()
        var_fp.wait()
        var_fpo.close()


def main():
    parser = ArgumentParser(description="Extract variant type and allele from a truth dataset")

    parser.add_argument('--vcf_fn', type=str, default="input.vcf", required=True,
                        help="Truth VCF file input, required")

    parser.add_argument('--var_fn', type=str, default="PIPE",
                        help="Truth variants output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help="Truth VCF file input, only used when vcf_fn is unified vcf. Marked truth variants not in unified as missing")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    OutputVariant(args)


if __name__ == "__main__":
    main()
