import sys
import shlex
from argparse import ArgumentParser
import shared.param_p as param
from shared.utils import subprocess_popen

def split_extend_bed(args):

    """
    Split bed file regions according to the contig name and extend bed region with no_of_positions =
    flankingBaseNum + 1 + flankingBaseNum, which allow samtools mpileup submodule to scan the flanking windows.
    """

    bed_fn = args.bed_fn
    output_fn = args.output_fn
    contig_name = args.ctgName
    region_start = args.ctgStart
    region_end = args.ctgEnd
    expand_region_size = args.expand_region_size
    if bed_fn is None:
        return
    output = []
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_fn)))
    pre_end, pre_start = -1, -1

    for row in unzip_process.stdout:

        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_name != None and ctg_name != contig_name:
            continue
        ctg_start, ctg_end = int(columns[1]), int(columns[2])
        if region_start and ctg_end < region_start:
            continue
        if region_end and ctg_start > region_end:
            break
        if pre_start == -1:
            pre_start = ctg_start - expand_region_size
            pre_end = ctg_end + expand_region_size
            continue
        if pre_end >= ctg_start - expand_region_size:
            pre_end = ctg_end + expand_region_size
            continue
        else:
            output.append(' '.join([contig_name, str(pre_start), str(pre_end)]))
            pre_start = ctg_start - expand_region_size
            pre_end = ctg_end + expand_region_size

    with open(output_fn, 'w') as output_file:
        output_file.write('\n'.join(output))

    unzip_process.stdout.close()
    unzip_process.wait()


def main():
    parser = ArgumentParser(description="Extend bed region for pileup calling")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Path to directory that stores small bins, default: %(default)s)"
                        )
    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Path of the output folder, default: %(default)s")

    parser.add_argument('--expand_region_size', type=int, default=param.no_of_positions,
                        help="Expand region size for each bed region, default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    split_extend_bed(args)


if __name__ == "__main__":
    main()