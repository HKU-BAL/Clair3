import os
import sys
from os.path import isfile, abspath
from sys import exit, stderr
from subprocess import check_output, PIPE, Popen
import argparse
import shlex
from subprocess import PIPE
from os.path import isfile, isdir
# A->A
# C->C
# G->G
# T or U->T
# R->A or G
# Y->C or T
# S->G or C
# W->A or T
# K->G or T
# M->A or C
# B->C or G or T
# D->A or G or T
# H->A or C or T
# V->A or C or G

def convert_iupac_to_n(string):
    if string == ".":
        return string
    output_str = []
    not_acgt_count = 0
    for s in string:
        if s.upper() not in "ACGTN,.":
            not_acgt_count += 1
            output_str.append('N')
        else:
            output_str.append(s)
    if not_acgt_count == 0:
        return string
    return ''.join(output_str)

IUPAC_base_to_ACGT_base_dict = dict(zip(
    "ACGTURYSWKMBDHVN",
    ("A", "C", "G", "T", "T", "A", "C", "C", "A", "G", "A", "C", "A", "A", "A", "A")
))

IUPAC_base_to_num_dict = dict(zip(
    "ACGTURYSWKMBDHVN",
    (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0)
))
BASIC_BASES = set("ACGTU")

WARNING = '\033[93m'
ERROR = '\033[91m'
ENDC = '\033[0m'

def log_error(log):
    return ERROR + log + ENDC

def log_warning(log):
    return WARNING + log + ENDC

def is_file_exists(file_name, suffix=""):
    if not isinstance(file_name, str) or not isinstance(suffix, str):
        return False
    return isfile(file_name + suffix)

def is_folder_exists(folder_name, suffix=""):
    if not isinstance(folder_name, str) or not isinstance(suffix, str):
        return False
    return isdir(folder_name + suffix)


def legal_range_from(param_name, x, min_num=None, max_num=None, exit_out_of_range=False):

    if min_num is not None and x < min_num and exit_out_of_range:
        exit(log_error("[ERROR] parameter --{}={} (minimum {}) out of range".format(param_name, x, min_num)))
    if max_num is not None and x > max_num and exit_out_of_range:
        exit(log_error("[ERROR] parameter --{}={} (maximum:{}) out of range".format(param_name, x, max_num)))
    return

def file_path_from(file_name, suffix="", exit_on_not_found=False, sep=""):
    if is_file_exists(file_name, suffix):
        return abspath(file_name + suffix)
    #allow fn.bam.bai->fn.bai fn.fa.fai->fn.fai
    elif sep != "" and len(sep) == 1:
        file_name_remove_suffix = sep.join(file_name.split(sep)[:-1])
        if is_file_exists(file_name_remove_suffix, suffix):
            return abspath(file_name_remove_suffix + suffix)
    if exit_on_not_found:
        exit(log_error("[ERROR] file %s not found" % (file_name + suffix)))
    return None

def folder_path_from(folder_name, create_not_found=True, exit_on_not_found=False):
    if is_folder_exists(folder_name):
        return abspath(folder_name)
    if exit_on_not_found:
        exit(log_error("[ERROR] folder %s not found" % (folder_name)))
    if create_not_found:
        if not os.path.exists(folder_name):
            os.makedirs(abspath(folder_name))
            print("[INFO] Create folder %s" % (folder_name), file=stderr)
            return abspath(folder_name)
    return None


def is_command_exists(command):
    if not isinstance(command, str):
        return False

    try:
        check_output("which %s" % (command), shell=True)
        return True
    except:
        return False


def executable_command_string_from(command_to_execute, exit_on_not_found=False):
    if is_command_exists(command_to_execute):
        return command_to_execute
    if exit_on_not_found:
        exit(log_error("[ERROR] %s executable not found" % (command_to_execute)))
    return None


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)

def str_none(v):
    if v is None:
        return None
    if v.upper() == "NONE":
        return None
    if isinstance(v, str):
        return v

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def region_from(ctg_name, ctg_start=None, ctg_end=None):
    """
    1-based region string [start, end]
    """
    if ctg_name is None:
        return ""
    if (ctg_start is None) != (ctg_end is None):
        return ""

    if ctg_start is None and ctg_end is None:
        return "{}".format(ctg_name)
    return "{}:{}-{}".format(ctg_name, ctg_start, ctg_end)

def reference_sequence_from(samtools_execute_command, fasta_file_path, regions):
    refernce_sequences = []
    region_value_for_faidx = " ".join(regions)

    samtools_faidx_process = subprocess_popen(
        shlex.split("{} faidx {} {}".format(samtools_execute_command, fasta_file_path, region_value_for_faidx))
    )
    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    # first line is reference name ">xxxx", need to be ignored
    reference_sequence = "".join(refernce_sequences[1:])

    # uppercase for masked sequences
    reference_sequence = reference_sequence.upper()

    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return reference_sequence

def vcf_candidates_from(vcf_fn, contig_name=None):

    known_variants_set =  set()
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))

    start_pos, end_pos = float('inf'), 0
    for row in unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split(maxsplit=3)
        ctg_name = columns[0]

        if contig_name and ctg_name != contig_name:
            continue
        center_pos = int(columns[1])
        known_variants_set.add(center_pos)
        start_pos = min(start_pos, center_pos)
        end_pos = max(center_pos, end_pos)

    known_variants_list = sorted(list(known_variants_set))
    return known_variants_list

def candidate_position_generator_from(
    candidate,
    flanking_base_num,
    begin_to_end
):
    for position in candidate:
        for i in range(position - (flanking_base_num + 1), position + (flanking_base_num + 1)):
            if i not in begin_to_end:
                begin_to_end[i] = [(position + (flanking_base_num + 1), position)]
            else:
                begin_to_end[i].append((position + (flanking_base_num + 1), position))
        yield position
    yield -1


def samtools_mpileup_generator_from(
    candidate,
    flanking_base_num,
    begin_to_end
):
    for position in candidate:
        for i in range(position - (flanking_base_num + 1), position + (flanking_base_num + 1)):
            if i not in begin_to_end:
                begin_to_end[i] = [(position + (flanking_base_num + 1), position)]
            else:
                begin_to_end[i].append((position + (flanking_base_num + 1), position))
        yield position
    yield -1

def samtools_view_process_from(
    ctg_name,
    ctg_start,
    ctg_end,
    samtools,
    bam_file_path
):
    have_start_and_end_position = ctg_start != None and ctg_end != None
    region_str = ("%s:%d-%d" % (ctg_name, ctg_start, ctg_end)) if have_start_and_end_position else ctg_name

    return subprocess_popen(
        shlex.split("%s view -F 2318 %s %s" % (samtools, bam_file_path, region_str))
    )

def get_header(reference_file_path=None, cmd_fn=None, sample_name="SAMPLE", version='1.1.2', gvcf=False, return_contig_length=False):
    from textwrap import dedent

    contig_length_dict = {}
    if reference_file_path is None or not os.path.exists(reference_file_path):
        ref_header_str = ""
    else:
        ref_header_str = "##reference={}".format(reference_file_path)

    cmdline_str = ""
    if cmd_fn is not None and os.path.exists(cmd_fn):
            cmd_line = open(cmd_fn).read().rstrip()

            if cmd_line is not None and len(cmd_line) > 0:
                cmdline_str = "##cmdline={}".format(cmd_line)

    if gvcf:
        header = dedent("""\
            ##fileformat=VCFv4.2
            ##source=Clair3
            ##clair3_version={}
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=LowQual,Description="Low quality variant">
            ##FILTER=<ID=RefCall,Description="Reference call">
            ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
            ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
            ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
            ##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
            ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">\n""".format(
            version))
    else:
        header = dedent("""\
            ##fileformat=VCFv4.2
            ##source=Clair3
            ##clair3_version={}
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=LowQual,Description="Low quality variant">
            ##FILTER=<ID=RefCall,Description="Reference call">
            ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
            ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
            ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">\n""".format(
            version))
    if ref_header_str != "":
        header_list = header.rstrip('\n').split('\n')
        insert_index = 3 if len(header_list) >= 3 else len(header_list) - 1
        header_list.insert(insert_index, ref_header_str)
        header = "\n".join(header_list) + '\n'

    if cmdline_str != "":
        header_list = header.rstrip('\n').split('\n')
        insert_index = 3 if len(header_list) >= 3 else len(header_list) - 1
        header_list.insert(insert_index, cmdline_str)
        header = "\n".join(header_list) + '\n'

    if reference_file_path is not None:
        reference_index_file_path = file_path_from(reference_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
        with open(reference_index_file_path, "r") as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")
                contig_name, contig_size = columns[0], columns[1]
                contig_length_dict[contig_name] = int(contig_size)
                header += "##contig=<ID=%s,length=%s>\n" % (contig_name, contig_size)

        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name)

    if return_contig_length:
        return header, contig_length_dict
    return header


