import os
import sys
import argparse
import shlex
from collections import defaultdict
from argparse import SUPPRESS
import shared.param_p as param
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import file_path_from, executable_command_string_from, folder_path_from, subprocess_popen

MIN_CHUNK_LENGTH = 200000
MAX_CHUNK_LENGTH = 20000000
major_contigs = {"chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]}.union(
    {str(a) for a in list(range(1, 23)) + ["X", "Y"]})
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def split_extend_bed(bed_fn, output_fn, contig_set):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_fn)))
    for row in unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if ctg_name not in contig_set:
            continue

        ctg_start, ctg_end = int(columns[1]), int(columns[2])
        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(ctg_start - expand_region_size), str(ctg_end + expand_region_size)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()


def CheckEnvs(args):
    basedir = os.path.dirname(__file__)
    bam_fn = file_path_from(args.bam_fn, exit_on_not_found=True)
    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    fai_fn = file_path_from(args.ref_fn + ".fai", exit_on_not_found=True)
    bai_fn = file_path_from(args.bam_fn + ".bai", exit_on_not_found=True)
    bed_fn = file_path_from(args.bed_fn)
    tree = bed_tree_from(bed_file_path=bed_fn)

    # create temp file folder
    output_fn_prefix = args.output_fn_prefix
    output_fn_prefix = folder_path_from(output_fn_prefix, create_not_found=True)
    log_path = folder_path_from(os.path.join(output_fn_prefix, 'log'), create_not_found=True)
    tmp_file_path = folder_path_from(os.path.join(output_fn_prefix, 'tmp'), create_not_found=True)
    split_bed_path = folder_path_from(os.path.join(tmp_file_path, 'split_beds'),
                                      create_not_found=True) if bed_fn else None
    pileup_vcf_path = folder_path_from(os.path.join(tmp_file_path, 'pileup_output'), create_not_found=True)
    merge_vcf_path = folder_path_from(os.path.join(tmp_file_path, 'merge_output'), create_not_found=True)
    phase_output_path = folder_path_from(os.path.join(tmp_file_path, 'phase_output'), create_not_found=True)
    gvcf_temp_output_path = folder_path_from(os.path.join(tmp_file_path, 'gvcf_tmp_output'), create_not_found=True)
    full_alignment_output_path = folder_path_from(os.path.join(tmp_file_path, 'full_alignment_output'),
                                                  create_not_found=True)
    phase_vcf_path = folder_path_from(os.path.join(phase_output_path, 'phase_vcf'), create_not_found=True)
    phase_bam_path = folder_path_from(os.path.join(phase_output_path, 'phase_bam'), create_not_found=True)
    candidate_bed_path = folder_path_from(os.path.join(full_alignment_output_path, 'candidate_bed'),
                                          create_not_found=True)

    is_include_all_contigs = args.includingAllContigs
    is_bed_file_provided = bed_fn is not None

    ctg_name_list = args.ctg_name
    is_ctg_name_list_provided = ctg_name_list is not None and ctg_name_list != "EMPTY"
    contig_set = set(ctg_name_list.split(',')) if is_ctg_name_list_provided else set()
    contig_set = contig_set.intersection(
        set(tree.keys())) if is_bed_file_provided and is_ctg_name_list_provided else contig_set
    contig_set = contig_set.union(
        set(tree.keys())) if is_bed_file_provided and not is_ctg_name_list_provided else contig_set

    # if each split region is too small(long) for given default chunk num, will increase(decrease) the total chunk num
    default_chunk_num = args.chunk_num
    DEFAULT_CHUNK_SIZE = args.chunk_size
    contig_length_list = []
    contig_chunk_num = {}
    with open(fai_fn, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")
            contig_name, contig_length = columns[0], int(columns[1])
            if not is_include_all_contigs and str(contig_name) not in major_contigs:
                continue

            if is_bed_file_provided and contig_name not in tree:
                continue
            if is_ctg_name_list_provided and contig_name not in contig_set:
                continue
            contig_set.add(contig_name)
            contig_length_list.append(contig_length)
            chunk_num = int(
                contig_length / float(DEFAULT_CHUNK_SIZE)) + 1 if contig_length % DEFAULT_CHUNK_SIZE else int(
                contig_length / float(DEFAULT_CHUNK_SIZE))
            contig_chunk_num[contig_name] = max(chunk_num, 1)


    if default_chunk_num > 0:
        min_chunk_length = min(contig_length_list) / float(default_chunk_num)
        max_chunk_length = max(contig_length_list) / float(default_chunk_num)
    if default_chunk_num > 0 and max_chunk_length > MAX_CHUNK_LENGTH:
        print(
        '[WARNING] Current maximum chunk size {} is larger than default maximum chunk size {}, You may set a larger chunk_num by setting --chunk_num=$ for better parallelism.'.format(
            min_chunk_length, MAX_CHUNK_LENGTH))

    elif default_chunk_num > 0 and min_chunk_length < MIN_CHUNK_LENGTH:
        print(
        '[WARNING] Current minimum chunk size {} is smaller than default minimum chunk size {}, You may set a smaller chunk_num by setting --chunk_num=$ .'.format(
            min_chunk_length, MIN_CHUNK_LENGTH))

    if default_chunk_num == 0 and max(contig_length_list) < DEFAULT_CHUNK_SIZE / 5:
        print(
        '[WARNING] Current maximum contig length {} is much smaller than default chunk size {}, You may set a smaller chunk size by setting --chunk_size=$ for better parallelism.'.format(
            max(contig_length_list), DEFAULT_CHUNK_SIZE))

    contigs_order = major_contigs_order + list(contig_set)

    sorted_contig_list = sorted(list(contig_set), key=lambda x: contigs_order.index(x))

    if not len(contig_set):
        exit("[ERROR] No contig name provide by --ctg_name or --bed_fn or provide contig name not found")
    else:
        print('[INFO] Call variant in contigs: {}'.format(' '.join(sorted_contig_list)))
        print('[INFO] Chunk number for each contig: {}'.format(
            ' '.join([str(contig_chunk_num[c]) for c in sorted_contig_list])))

    if is_bed_file_provided:
        split_extend_bed(bed_fn=bed_fn, output_fn=split_bed_path, contig_set=contig_set)

    contig_name_list = os.path.join(tmp_file_path, 'CONTIGS')
    chunk_list = os.path.join(tmp_file_path, 'CHUNK_LIST')

    with open(contig_name_list, 'w') as output_file:
        output_file.write('\n'.join(sorted_contig_list))

    with open(chunk_list, 'w') as output_file:
        for contig_name in sorted_contig_list:
            chunk_num = contig_chunk_num[contig_name]
            for chunk_id in range(1, chunk_num + 1):
                output_file.write(contig_name + ' ' + str(chunk_id) + ' ' + str(chunk_num) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Check the envorinment variables and handle contig and bed if need")

    parser.add_argument('--bam_fn', type=str, default=None,
                        help="BAM file input, default: %(default)s")

    parser.add_argument('--output_fn_prefix', type=str, default=None,
                        help="Path to the output folder")

    parser.add_argument('--ctg_name', type=str, default='EMPTY',
                        help="The name of contigs to be processed")

    parser.add_argument('--bed_fn', type=str, nargs='?', action="store", default=None,
                        help="Call variant only in these regions, works in intersection with ctgName, ctgStart and ctgEnd, optional, default: as defined by ctgName, ctgStart and ctgEnd")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--chunk_size', type=int, default=5000000,
                        help="Total chunk number for parallel execution. Each chunk refer to a smaller reference regions.")

    parser.add_argument('--includingAllContigs', action='store_true',
                        help="Call variants on all contigs, default: chr{1..22,X,Y,M,MT} and {1..22,X,Y,MT}")

    # options for internal process control
    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=0,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    if not args.includingAllContigs and args.ctg_name == 'EMPTY':
        print("[INFO] --includingAllContigs not enabled, use chr{1..22,X,Y,M,MT} and {1..22,X,Y,MT} by default")
    else:
        print("[INFO] --includingAllContigs enabled")

    CheckEnvs(args)


if __name__ == "__main__":
    main()
