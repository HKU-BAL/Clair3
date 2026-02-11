import os
import sys
import shlex
import subprocess
import argparse
import concurrent.futures

from collections import Counter
from sys import stderr
from subprocess import PIPE, run
from argparse import ArgumentParser
from collections import defaultdict
from subprocess import Popen

def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_sv_qual(input_variant_dict, tree, contig, ctg_name, region_start):

    if not tree or (contig is None) or (contig not in tree):
        return None

    interval_tree = tree[contig]
    pos = int(list(interval_tree.at(region_start))[0].begin)
    key = pos if ctg_name is not None else (contig, pos)
    qual = round(float(input_variant_dict[key].qual),2) if key in input_variant_dict else None

    return qual



def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


class Position(object):
    def __init__(self,
                 ctg_name=None,
                 genotype1=None,
                 genotype2=None,
                 pos=None,
                 ref_base=None,
                 alt_base=None,
                 filter=None,
                 depth=None,
                 af=None,
                 qual=None,
                 sv_del_end=None,
                 row_str=None):
        self.ctg_name = ctg_name
        self.pos = pos
        self.reference_bases = ref_base

        self.alternate_bases = [alt_base] if ',' not in alt_base else alt_base.split(',')
        self.genotype1 = genotype1
        self.genotype2 = genotype2
        self.genotype = [genotype1, genotype2]
        self.depth = depth
        self.af = af
        self.qual = qual
        self.row_str = row_str
        self.filter = filter
        self.sv_del_end = sv_del_end


class VcfWriter(object):
    def __init__(self, vcf_fn, ctg_name=None, sample_name="SAMPLE", show_ref_calls=False):
        self.vcf_fn = vcf_fn
        self.show_ref_calls = show_ref_calls

        # make directory if not exist
        vcf_folder = os.path.dirname(self.vcf_fn)
        if not os.path.exists(vcf_folder):
            print("[INFO] Output VCF folder {} not found, create it".format(vcf_folder))
            return_code = run("mkdir -p {}".format(vcf_folder), shell=True)

        self.vcf_writer = open(self.vcf_fn, 'w')
        self.ctg_name = ctg_name
        self.sample_name = sample_name

    def close(self):
        try:
            self.vcf_writer.close()
        except:
            pass

class VcfReader(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 show_ref=True,
                 direct_open=False,
                 keep_row_str=False,
                 save_header=False,
                 filter_tag=None,
                 sv_input=False,
                 sv_alt_tag=None,
                 ):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.variant_dict = defaultdict(Position)
        self.show_ref = show_ref
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.header = ""
        self.filter_tag = filter_tag
        self.save_header = save_header
        self.sv_input = sv_input
        self.sv_alt_tag = sv_alt_tag

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
            sv_del_end = None
            FILTER = None
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue

            if self.filter_tag is not None:
                FILTER = columns[6]
                filter_list = self.filter_tag.split(',')
                if sum([1 if filter == FILTER else 0 for filter in filter_list]) == 0:
                    continue

            reference, alternate, next_to_last_column, last_column = columns[3], columns[4], columns[-2], columns[-1]

            if self.sv_input:
                if self.sv_alt_tag is not None:
                    if alternate != self.sv_alt_tag:
                        continue
                try:
                    INFO = columns[7].split(';')
                    for info in INFO:
                        if info.startswith('END='):
                            sv_del_end = int(info[4:])
                            if sv_del_end < position:
                                continue
                except:
                    continue

            qual = columns[5] if len(columns) > 5 else None

            genotype = last_column.split(":")[0].replace("/", "|").replace(".", "0").split("|")
            try:
                genotype_1, genotype_2 = genotype
                if int(genotype_1) > int(genotype_2):
                    genotype_1, genotype_2 = genotype_2, genotype_1
            except:
                genotype_1 = -1
                genotype_2 = -1

            # remove * to guarentee vcf match
            if '*' in alternate:
                alternate = alternate.split(',')
                if int(genotype_1) + int(genotype_2) != 3 or len(alternate) != 2:
                    continue

                alternate = ''.join([alt_base for alt_base in alternate if alt_base != '*'])

            try:
                af_idx = next_to_last_column.split(':').index('AF')
                af = float(last_column.split(':')[af_idx])
            except:
                af = None

            if genotype_1 == "0" and genotype_2 == "0" and not self.show_ref:
                continue

            row_str = row if self.keep_row_str else False
            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              ref_base=reference,
                                              alt_base=alternate,
                                              genotype1=int(genotype_1),
                                              genotype2=int(genotype_2),
                                              qual=qual,
                                              af=af,
                                              filter=FILTER,
                                              sv_del_end=sv_del_end,
                                              row_str=row_str)


def get_base_list(columns):
    pileup_bases = columns[4]

    base_idx = 0
    base_list = []
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx].upper()
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]  # add indel seq
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':  # start of read
            base_idx += 1
        # skip $, the end of read
        base_idx += 1

    upper_base_counter = Counter([''.join(item).upper() for item in base_list])
    return upper_base_counter, base_list


def update_header(header):
    header = header.rstrip().split("\n")
    header.insert(-1,
                  '##FORMAT=<ID=PEAD,Number=1,Type=Integer,Description="Allelic depths for the alt alleles after removing same pair-end reads in the order listed"')
    return '\n'.join(header) + '\n'


def extract_base(POS):
    pos = POS.pos
    bam_fn = args.bam_fn
    samtools = args.samtools
    ctg_name = args.ctg_name if args.ctg_name is not None else POS.ctg_name

    min_mq = args.min_mq
    min_bq = args.min_bq

    ref_base = POS.reference_bases
    alt_base = POS.alternate_bases[0]
    ctg_range = "{}:{}-{}".format(ctg_name, pos, pos)
    samtools_command = "{} mpileup {} --min-MQ {} --min-BQ {} --excl-flags 2316 -r {} --output-QNAME".format(samtools, bam_fn, min_mq,
                                                                                              min_bq, ctg_range)

    output = subprocess.run(samtools_command, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    output = output.stdout.rstrip()

    key = (ctg_name, int(pos))
    for row in output.rstrip().split('\n'):
        columns = row.split('\t')
        ctg_name, p = columns[0], int(columns[1])
        if p != POS.pos:
            continue

        base_counter, base_list = get_base_list(columns)
        read_name_list = columns[6].split(',')
        match_read_read_name = []
        match_index = []
        for idx, (base, rn) in enumerate(zip(base_list, read_name_list)):
            if len(ref_base) == 1 and len(alt_base) == 1:
                if base[0].upper() == alt_base and base[1] == "":
                    match_index.append(idx)
                    match_read_read_name.append(rn)
            elif len(ref_base) == 1 and len(alt_base) > 1:
                if ''.join(base).replace('+', '').upper() == alt_base and '+' in ''.join(base):
                    match_index.append(idx)
                    match_read_read_name.append(rn)
            elif len(ref_base) > 1 and len(alt_base) == 1:
                if base[0].upper() == ref_base[0] and '-' in base[1] and len(base[1]) == (len(ref_base)):
                    match_index.append(idx)
                    match_read_read_name.append(rn)

    return key, len(set(match_read_read_name))


def Run(args):
    ctg_name = args.ctg_name
    threads = args.threads
    clair3_vcf_input = args.clair3_vcf_input
    filter_tag = args.filter_tag


    input_vcf_reader = VcfReader(vcf_fn=clair3_vcf_input,
                                 ctg_name=ctg_name,
                                 show_ref=False,
                                 keep_row_str=True,
                                 filter_tag=filter_tag,
                                 save_header=True)

    input_vcf_reader.read_vcf()
    input_variant_dict = input_vcf_reader.variant_dict

    vcf_output = args.vcf_output
    if vcf_output.endswith('.gz'):
        vcf_output = vcf_output[:-3]
    vcf_writer = VcfWriter(vcf_fn=vcf_output, ctg_name=ctg_name, show_ref_calls=False)


    result_dict = defaultdict(int)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as exec:
        for result in exec.map(extract_base, list(input_variant_dict.values())):
           result_dict[result[0]] = result[1]

    vcf_writer.vcf_writer.write(update_header(input_vcf_reader.header))
    for key in input_variant_dict.keys():
        pos = key if ctg_name is not None else key[1]
        contig = ctg_name if ctg_name is not None else key[0]
        row_str = input_variant_dict[key].row_str.rstrip()
        key = (contig, pos)
        if key in result_dict:
            columns = row_str.rstrip().split('\t')
            columns[8] += ':' + "PEAD"
            columns[9] += ':' + str(result_dict[key])
            row_str = '\t'.join(columns)

        vcf_writer.vcf_writer.write(row_str + '\n')

    vcf_writer.close()

    compress_index_vcf(vcf_output)


def main():
    parser = ArgumentParser(description="Add pair-end allele depth to Clair3 VCF file output")

    parser.add_argument('--bam_fn', type=str, default=None,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--clair3_vcf_input', type=str, default=None,
                        help="Clair3 VCF input")

    parser.add_argument('--vcf_output', type=str, default=None,
                        help="VCF output path")

    parser.add_argument('--threads', type=int, default=8,
                        help="Max #threads to be used")

    parser.add_argument('--filter_tag', type=str, default="PASS",
                        help="Filter tag for the input VCF")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--min_mq', type=int, default=5,
                        help="DEBUG: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=0,
                        help="DEBUG: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")


    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)
    global args

    args = parser.parse_args()

    Run(args)


if __name__ == "__main__":
    main()


