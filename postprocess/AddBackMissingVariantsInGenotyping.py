import os
import subprocess
import shlex
import argparse

from sys import stderr
from subprocess import PIPE, run
from argparse import ArgumentParser
from collections import defaultdict
from subprocess import Popen

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'ture', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'flase', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


class Position(object):
    def __init__(self,
                 ctg_name=None,
                 pos=None,
                 row_str=None):
        self.ctg_name = ctg_name
        self.pos = pos
        self.row_str = row_str


class VcfReader(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 direct_open=False,
                 keep_row_str=False,
                 save_header=False,
                 ):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.variant_dict = defaultdict(Position)
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.header = ""
        self.save_header = save_header

    def read_vcf(self):

        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        if self.direct_open:
            vcf_fp = open(self.vcf_fn)
            vcf_fo = vcf_fp
        else:
            vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
            vcf_fo = vcf_fp.stdout
        for row in vcf_fo:
            columns = row.strip().split()
            if columns[0][0] == "#":
                if self.save_header:
                    self.header += row
                continue

            # position in vcf is 1-based
            chromosome, position = columns[0], int(columns[1])
            key = (chromosome, position) if self.ctg_name is None else position
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue

            row_str = row if self.keep_row_str else False
            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              row_str=row_str)


def genotype_vcf(args):
    vcf_fn = args.vcf_fn
    clair3_input_vcf_fn = args.clair3_input_vcf_fn
    output_fn = args.output_fn
    switch_genotype = args.switch_genotype

    output = open(output_fn, 'w')
    vcf_reader = VcfReader(vcf_fn=vcf_fn,
                           ctg_name=None,
                           keep_row_str=True,
                           save_header=True)

    vcf_reader.read_vcf()
    variant_dict = vcf_reader.variant_dict

    clair3_vcf_reader = VcfReader(vcf_fn=clair3_input_vcf_fn,
                                  ctg_name=None,
                                  keep_row_str=True,
                                  save_header=True)

    clair3_vcf_reader.read_vcf()
    clair3_variant_dict = clair3_vcf_reader.variant_dict

    all_contigs_list = list(set([k[0] for k in variant_dict]))

    contigs_order = major_contigs_order + all_contigs_list

    contigs_order_list = sorted(all_contigs_list, key=lambda x: contigs_order.index(x))

    output.write(clair3_vcf_reader.header)



    count = 0
    contig_dict = defaultdict(list)
    for k, v in variant_dict.items():
        ctg, pos = k
        if k not in clair3_variant_dict:
            row_str = variant_dict[k].row_str
            count += 1
            if switch_genotype:
                columns = row_str.rstrip().split('\t')
                if len(columns) < 10:
                    columns += ['.'] * (10 - len(columns))
                columns[3] = columns[3][0] if len(columns[3]) > 0 else '.' # Only keep the reference base for REF
                columns[4] = '.' # ALT to 0
                columns[5] = "." # QUAL to .
                columns[6] = "." # FILTER to .
                columns[7] = "." # INFO to .
                columns[8] = "GT" # keep GT tag only
                columns[9] = './.'
                row_str = '\t'.join(columns) + '\n'
        else:
            row_str = clair3_variant_dict[k].row_str

        contig_dict[ctg].append((int(pos), row_str))

    for contig in contigs_order_list:
        row_list = [item[1] for item in sorted(contig_dict[contig], key=lambda x: x[0])]
        output.write(''.join(row_list))

    output.close()

    print("[INFO] Total variants for genotyping: {}, total Clair3 variant output: {}, added {} variants into output VCF"\
          .format(len(variant_dict), len(clair3_variant_dict), count))
    compress_index_vcf(output_fn)


def main():
    parser = ArgumentParser(description="Genotype VCF in postprocessing")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input for genotyping")

    parser.add_argument('--clair3_input_vcf_fn', type=str, default=None,
                        help="Clair3 input vcf")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Output vcf file name")

    parser.add_argument('--switch_genotype', type=str2bool, default=True,
                        help="Switch missed variant genotype to ./.")

    args = parser.parse_args()

    genotype_vcf(args)


if __name__ == "__main__":
    main()
