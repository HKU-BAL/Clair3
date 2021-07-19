import os
import sys
import argparse
import shlex
import subprocess

from collections import defaultdict
from argparse import SUPPRESS
from distutils.version import LooseVersion

import shared.param_p as param
from shared.interval_tree import bed_tree_from
from shared.utils import file_path_from, folder_path_from, subprocess_popen, str2bool, \
    legal_range_from, log_error, log_warning

MIN_CHUNK_LENGTH = 200000
MAX_CHUNK_LENGTH = 20000000
major_contigs = {"chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]}.union(
    {str(a) for a in list(range(1, 23)) + ["X", "Y"]})
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

required_tool_version = {
    'python': LooseVersion('3.6.10'),
    'pypy': LooseVersion('3.6'),
    'samtools': LooseVersion('1.10'),
    'whatshap': LooseVersion('1.0'),
    'parallel': LooseVersion('20191122'),
}

def check_version(tool, pos=None, is_pypy=False):
    try:
        if is_pypy:
            proc = subprocess.run("{} -c 'import sys; print (sys.version)'".format(tool), stdout=subprocess.PIPE,
                                  shell=True)
        else:
            proc = subprocess.run([tool, "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        first_line = proc.stdout.decode().split("\n", 1)[0]
        version = first_line.split()[pos]
        version = LooseVersion(version)
    except Exception:
        return None

    return version


def check_python_path():
    python_path = subprocess.run("which python", stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
    sys.exit(log_error("[ERROR] Current python execution path: {}".format(python_path)))


def check_tools_version(tool_version, required_tool_version):
    for tool, version in tool_version.items():
        required_version = required_tool_version[tool]
        if version is None:
            print(log_error("[ERROR] {} not found, please check you are in clair3 virtual environment".format(tool)))
            check_python_path()
        elif version < required_version:
            print(log_error("[ERROR] Tool version not match, please check you are in clair3 virtual environment"))
            print(' '.join([str(item).ljust(10) for item in ["Tool", "Version", "Required"]]))
            error_info = ' '.join([str(item).ljust(10) for item in [tool, version, '>=' + str(required_version)]])
            print(error_info)
            check_python_path()
    return


def check_contig_in_bam(bam_fn, sorted_contig_list, samtools):
    bai_process = subprocess_popen(shlex.split("{} idxstats {}".format(samtools, bam_fn)))
    contig_with_read_support_set = set()
    for row_id, row in enumerate(bai_process.stdout):
        row = row.split('\t')
        if len(row) != 4:
            continue
        contig_name, contig_length, mapped_reads, unmapped_reads = row
        if contig_name not in sorted_contig_list:
            continue
        if int(mapped_reads) > 0:
            contig_with_read_support_set.add(contig_name)
    for contig_name in sorted_contig_list:
        if contig_name not in contig_with_read_support_set:
            print(log_warning(
                "[WARNING] Contig name {} provided but no mapped reads in BAM, skip!".format(contig_name)))
    filtered_sorted_contig_list = [item for item in sorted_contig_list if item in contig_with_read_support_set]

    found_contig = True
    if len(filtered_sorted_contig_list) == 0:
        found_contig = False
        print(log_warning(
            "[WARNING] No mapped reads support in BAM for provided contigs set {}".format(
                ' '.join(sorted_contig_list))))
    return filtered_sorted_contig_list, found_contig


def split_extend_vcf(vcf_fn, output_fn):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))

    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split(maxsplit=3)
        ctg_name = columns[0]

        center_pos = int(columns[1])
        ctg_start, ctg_end = center_pos - 1, center_pos
        if ctg_start < 0:
            sys.exit(
                log_error("[ERROR] Invalid VCF input in {}-th row {} {} {}".format(row_id + 1, ctg_name, center_pos)))
        if ctg_start - expand_region_size < 0:
            continue
        expand_ctg_start = ctg_start - expand_region_size
        expand_ctg_end = ctg_end + expand_region_size

        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(expand_ctg_start), str(expand_ctg_end)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()

    know_vcf_contig_set = set(list(output_ctg_dict.keys()))

    return know_vcf_contig_set


def split_extend_bed(bed_fn, output_fn, contig_set=None):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_fn)))
    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_set and ctg_name not in contig_set:
            continue

        ctg_start, ctg_end = int(columns[1]), int(columns[2])

        if ctg_end < ctg_start or ctg_start < 0 or ctg_end < 0:
            sys.exit(log_error(
                "[ERROR] Invalid BED input in {}-th row {} {} {}".format(row_id + 1, ctg_name, ctg_start, ctg_end)))
        expand_ctg_start = max(0, ctg_start - expand_region_size)
        expand_ctg_end = max(0, ctg_end + expand_region_size)
        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(expand_ctg_start), str(expand_ctg_end)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()


def output_header(output_fn, reference_file_path, sample_name='SAMPLE'):
    output_file = open(output_fn, "w")
    from textwrap import dedent
    output_file.write(dedent("""\
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FILTER=<ID=LowQual,Description="Low quality variant">
        ##FILTER=<ID=RefCall,Description="Reference call">
        ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
        ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
        ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">"""
                             ) + '\n')

    if reference_file_path is not None:
        reference_index_file_path = file_path_from(reference_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
        with open(reference_index_file_path, "r") as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")
                contig_name, contig_size = columns[0], columns[1]
                output_file.write(("##contig=<ID=%s,length=%s>" % (contig_name, contig_size) + '\n'))

    output_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name))
    output_file.close()

def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

def CheckEnvs(args):
    basedir = os.path.dirname(__file__)
    bam_fn = file_path_from(args.bam_fn, exit_on_not_found=True)
    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    fai_fn = file_path_from(args.ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')
    bai_fn = file_path_from(args.bam_fn, suffix=".bai", exit_on_not_found=True, sep='.')
    bed_fn = file_path_from(args.bed_fn)
    vcf_fn = file_path_from(args.vcf_fn)
    tree = bed_tree_from(bed_file_path=bed_fn)

    # create temp file folder
    output_fn_prefix = args.output_fn_prefix
    output_fn_prefix = folder_path_from(output_fn_prefix, create_not_found=True)
    log_path = folder_path_from(os.path.join(output_fn_prefix, 'log'), create_not_found=True)
    tmp_file_path = folder_path_from(os.path.join(output_fn_prefix, 'tmp'), create_not_found=True)
    split_bed_path = folder_path_from(os.path.join(tmp_file_path, 'split_beds'),
                                      create_not_found=True) if bed_fn or vcf_fn else None
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

    # environment parameters
    pypy = args.pypy
    samtools = args.samtools
    whatshap = args.whatshap
    parallel = args.parallel
    qual = args.qual
    var_pct_full = args.var_pct_full
    ref_pct_full = args.ref_pct_full
    snp_min_af = args.snp_min_af
    indel_min_af = args.indel_min_af
    sample_name = args.sampleName
    contig_name_list = os.path.join(tmp_file_path, 'CONTIGS')
    chunk_list = os.path.join(tmp_file_path, 'CHUNK_LIST')

    legal_range_from(param_name="qual", x=qual, min_num=0, exit_out_of_range=True)
    legal_range_from(param_name="var_pct_full", x=var_pct_full, min_num=0, max_num=1, exit_out_of_range=True)
    legal_range_from(param_name="ref_pct_full", x=ref_pct_full, min_num=0, max_num=1, exit_out_of_range=True)
    legal_range_from(param_name="snp_min_af", x=snp_min_af, min_num=0, max_num=1, exit_out_of_range=True)
    legal_range_from(param_name="indel_min_af", x=indel_min_af, min_num=0, max_num=1, exit_out_of_range=True)
    if ref_pct_full > 0.3:
        print(log_warning(
            "[WARNING] For efficiency, we use a maximum 30% reference candidates for full-alignment calling"))
    tool_version = {
        'python': LooseVersion(sys.version.split()[0]),
        'pypy': check_version(tool=pypy, pos=0, is_pypy=True),
        'samtools': check_version(tool=samtools, pos=1),
        'whatshap': check_version(tool=whatshap, pos=1),
        'parallel': check_version(tool=parallel, pos=2),
    }
    check_tools_version(tool_version, required_tool_version)

    is_include_all_contigs = args.include_all_ctgs
    is_bed_file_provided = bed_fn is not None
    is_known_vcf_file_provided = vcf_fn is not None

    if is_known_vcf_file_provided and is_bed_file_provided:
        sys.exit(log_error("[ERROR] Please provide either --vcf_fn or --bed_fn only"))

    if is_known_vcf_file_provided:
        know_vcf_contig_set = split_extend_vcf(vcf_fn=vcf_fn, output_fn=split_bed_path)

    ctg_name_list = args.ctg_name
    is_ctg_name_list_provided = ctg_name_list is not None and ctg_name_list != "EMPTY"
    contig_set = set(ctg_name_list.split(',')) if is_ctg_name_list_provided else set()

    if is_ctg_name_list_provided and is_bed_file_provided:
        print(log_warning("[WARNING] both --ctg_name and --bed_fn provided, will only proceed contigs in intersection"))

    if is_ctg_name_list_provided and is_known_vcf_file_provided:
        print(log_warning("[WARNING] both --ctg_name and --vcf_fn provided, will only proceed contigs in intersection"))

    if is_ctg_name_list_provided:

        contig_set = contig_set.intersection(
            set(tree.keys())) if is_bed_file_provided else contig_set

        contig_set = contig_set.intersection(
            know_vcf_contig_set) if is_known_vcf_file_provided else contig_set
    else:
        contig_set = contig_set.union(
            set(tree.keys())) if is_bed_file_provided else contig_set

        contig_set = contig_set.union(
            know_vcf_contig_set) if is_known_vcf_file_provided else contig_set

    # if each split region is too small(long) for given default chunk num, will increase(decrease) the total chunk num
    default_chunk_num = args.chunk_num
    DEFAULT_CHUNK_SIZE = args.chunk_size
    contig_length_list = []
    contig_chunk_num = {}

    threads = args.threads
    sched_getaffinity_list = list(os.sched_getaffinity(0))
    numCpus = len(sched_getaffinity_list)

    if threads > numCpus:
        print(log_warning(
            '[WARNING] Current maximum threads {} is larger than support cpu count {}, You may set a smaller parallel threads by setting --threads=$ for better parallelism.'.format(
                threads, numCpus)))

    ## for better parallelism for create tensor and call variants, we over commit the overall threads/4 for 3 times, which is 0.75 * overall threads.
    threads_over_commit = max(4, int(threads * 0.75))

    with open(fai_fn, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")
            contig_name, contig_length = columns[0], int(columns[1])
            if not is_include_all_contigs and (
            not (is_bed_file_provided or is_ctg_name_list_provided or is_known_vcf_file_provided)) and str(
                    contig_name) not in major_contigs:
                continue

            if is_bed_file_provided and contig_name not in tree:
                continue
            if is_ctg_name_list_provided and contig_name not in contig_set:
                continue
            if is_known_vcf_file_provided and contig_name not in contig_set:
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

    contigs_order = major_contigs_order + list(contig_set)

    sorted_contig_list = sorted(list(contig_set), key=lambda x: contigs_order.index(x))

    found_contig = True
    if not len(contig_set):
        if is_bed_file_provided:
            all_contig_in_bed = ' '.join(list(tree.keys()))
            print(log_warning("[WARNING] No contig intersection found by --bed_fn, contigs in BED {}: {}".format(bed_fn, all_contig_in_bed)))
        if is_known_vcf_file_provided:
            all_contig_in_vcf = ' '.join(list(know_vcf_contig_set))
            print(log_warning("[WARNING] No contig intersection found by --vcf_fn, contigs in VCF {}: {}".format(vcf_fn, all_contig_in_vcf)))
        if is_ctg_name_list_provided:
            all_contig_in_ctg_name = ' '.join(ctg_name_list.split(','))
            print(log_warning("[WARNING] No contig intersection found by --ctg_name, contigs in contigs list: {}".format(all_contig_in_ctg_name)))
        found_contig = False
    else:
        for c in sorted_contig_list:
            if c not in contig_chunk_num:
                print(log_warning(("[WARNING] Contig {} given but not found in reference fai file".format(c))))

        # check contig in bam have support reads
        sorted_contig_list, found_contig = check_contig_in_bam(bam_fn=bam_fn, sorted_contig_list=sorted_contig_list,
                                                               samtools=samtools)

    if not found_contig:
        # output header only to merge_output.vcf.gz
        output_fn = os.path.join(output_fn_prefix, "merge_output.vcf")
        output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
        compress_index_vcf(output_fn)
        print(log_warning(
            ("[WARNING] No contig intersection found, output header only in {}").format(output_fn + ".gz")))
        with open(contig_name_list, 'w') as output_file:
            output_file.write("")
        return

    print('[INFO] Call variant in contigs: {}'.format(' '.join(sorted_contig_list)))
    print('[INFO] Chunk number for each contig: {}'.format(
        ' '.join([str(contig_chunk_num[c]) for c in sorted_contig_list])))

    if default_chunk_num > 0 and max_chunk_length > MAX_CHUNK_LENGTH:
        print(log_warning(
            '[WARNING] Current maximum chunk size {} is larger than default maximum chunk size {}, You may set a larger chunk_num by setting --chunk_num=$ for better parallelism.'.format(
                min_chunk_length, MAX_CHUNK_LENGTH)))

    elif default_chunk_num > 0 and min_chunk_length < MIN_CHUNK_LENGTH:
        print(log_warning(
            '[WARNING] Current minimum chunk size {} is smaller than default minimum chunk size {}, You may set a smaller chunk_num by setting --chunk_num=$.'.format(
                min_chunk_length, MIN_CHUNK_LENGTH)))

    if default_chunk_num == 0 and max(contig_length_list) < DEFAULT_CHUNK_SIZE / 5:
        print(log_warning(
            '[WARNING] Current maximum contig length {} is much smaller than default chunk size {}, You may set a smaller chunk size by setting --chunk_size=$ for better parallelism.'.format(
                max(contig_length_list), DEFAULT_CHUNK_SIZE)))

    if is_bed_file_provided:
        split_extend_bed(bed_fn=bed_fn, output_fn=split_bed_path, contig_set=contig_set)

    with open(contig_name_list, 'w') as output_file:
        output_file.write('\n'.join(sorted_contig_list))

    with open(chunk_list, 'w') as output_file:
        for contig_name in sorted_contig_list:
            chunk_num = contig_chunk_num[contig_name]
            for chunk_id in range(1, chunk_num + 1):
                output_file.write(contig_name + ' ' + str(chunk_id) + ' ' + str(chunk_num) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Check the environment and the validity of the input variables, preprocess the BED input if necessary")

    parser.add_argument('--bam_fn', type=str, default=None,
                        help="BAM file input, default: %(default)s")

    parser.add_argument('--output_fn_prefix', type=str, default=None,
                        help="Path to the output folder")

    parser.add_argument('--ctg_name', type=str, default='EMPTY',
                        help="The name of sequence to be processed, separated by comma")

    parser.add_argument('--bed_fn', type=str, nargs='?', action="store", default=None,
                        help="Call variant only in these regions. Will take an intersection if --ctg_name is set")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--chunk_size', type=int, default=5000000,
                        help="The size of each chuck for parallel processing, default: 5Mbp")

    parser.add_argument('--include_all_ctgs', type=str2bool, default=False,
                        help="Call variants on all contigs, default: chr{1..22,X,Y,M,MT} and {1..22,X,Y,MT}")

    parser.add_argument('--threads', type=int, default=16,
                        help="Max #threads to be used. The full genome will be divided into small chucks for parallel processing")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    parser.add_argument('--pypy', type=str, default="pypy3",
                        help="Path to the 'pypy', pypy3 version >= 3.6 is required, default: %(default)s")

    parser.add_argument('--python', type=str, default="python3",
                        help="Path to the 'python3', default: %(default)s")

    parser.add_argument('--parallel', type=str, default="parallel",
                        help="Path to the 'parallel', default: %(default)s")

    parser.add_argument('--whatshap', type=str, default="whatshap",
                        help="Path to the 'whatshap', default: %(default)s")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--qual', type=int, default=None,
                        help="If set, variants with >=$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--var_pct_full', type=float, default=0.3,
                        help="Default variant call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ref_pct_full', type=float, default=0.3,
                        help="Default reference call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--snp_min_af', type=float, default=0.08,
                        help="Minimum SNP allele frequency for a site to be considered as a candidate site, default: %(default)f")

    parser.add_argument('--indel_min_af', type=float, default=0.08,
                        help="Minimum Indel allele frequency for a site to be considered as a candidate site, default: %(default)f")

    # options for internal process control
    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=0,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    if not args.include_all_ctgs and args.ctg_name == 'EMPTY':
        print("[INFO] --include_all_ctgs not enabled, use chr{1..22,X,Y} and {1..22,X,Y} by default")
    elif args.include_all_ctgs:
        print("[INFO] --include_all_ctgs enabled")
        print(log_warning("[WARNING] Please enable --no_phasing_for_fa if calling variant in non-diploid organisms"))

    CheckEnvs(args)


if __name__ == "__main__":
    main()
