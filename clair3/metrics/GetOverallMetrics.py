import os
import sys
import shlex
import logging

from sys import stderr
from subprocess import Popen
from argparse import ArgumentParser
from subprocess import PIPE

logging.basicConfig(format='%(message)s', level=logging.INFO)


def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def metrics(query_fp, query_tp, truth_fn, truth_tp):
    # https://github.com/Illumina/hap.py/blob/master/doc/happy.md
    precision = query_tp / (query_tp + query_fp)
    recall = truth_tp / (truth_tp + truth_fn)
    f1_score = 2 * precision * recall / (precision + recall)
    return round(precision, 6), round(recall, 6), round(f1_score, 6)


def Cal(args):
    happy_vcf_fn = args.happy_vcf_fn
    contig_name = args.ctgName
    output_fn = args.output_fn

    if output_fn:
        output_file = open(output_fn, 'w')
    else:
        output_file = None
    happy_vcf_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (happy_vcf_fn)))

    truth_all_tp, query_all_tp, query_all_fp, truth_all_fn = 0, 0, 0, 0
    truth_snp_tp, query_snp_tp, query_snp_fp, truth_snp_fn = 0, 0, 0, 0
    truth_indel_tp, query_indel_tp, query_indel_fp, truth_indel_fn = 0, 0, 0, 0
    truth_ins_tp, query_ins_tp, query_ins_fp, truth_ins_fn = 0, 0, 0, 0
    truth_del_tp, query_del_tp, query_del_fp, truth_del_fn = 0, 0, 0, 0

    for row in happy_vcf_unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()

        ctg_name, pos = columns[0], int(columns[1])
        if contig_name is not None and ctg_name != contig_name:
            continue

        FORMAT, TRUTH, QUERY = columns[8], columns[9], columns[10]
        FORMAT = FORMAT.split(':')
        TRUTH = TRUTH.split(':')
        QUERY = QUERY.split(':')

        ft_dict = dict(zip(FORMAT, TRUTH))
        fq_dict = dict(zip(FORMAT, QUERY))

        # hap.py vcf header
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision for call (TP/FP/FN/N)">
        ##FORMAT=<ID=BK,Number=1,Type=String,Description="Sub-type for decision (match/mismatch type). (Loose match distance is 30)">
        ##FORMAT=<ID=BI,Number=1,Type=String,Description="Additional comparison information">
        ##FORMAT=<ID=QQ,Number=1,Type=Float,Description="Variant quality for ROC creation">
        ##INFO=<ID=Regions,Number=.,Type=String,Description="Tags for regions.">
        ##FORMAT=<ID=BVT,Number=1,Type=String,Description="High-level variant type (SNP|INDEL).">
        ##FORMAT=<ID=BLT,Number=1,Type=String,Description="High-level location type (het|homref|hetalt|homalt|nocall).">

        t_BD = ft_dict['BD'] if 'BD' in ft_dict else None
        t_BI = ft_dict['BI'] if 'BI' in ft_dict else None
        t_BVT = ft_dict['BVT'] if 'BVT' in ft_dict else None
        q_BD = fq_dict['BD'] if 'BD' in fq_dict else None
        q_BI = fq_dict['BI'] if 'BI' in fq_dict else None
        q_BVT = fq_dict['BVT'] if 'BVT' in fq_dict else None
        if not t_BD or not t_BI or not t_BVT or not q_BD or not q_BI or not q_BVT:
            sys.exit("[ERROR] Happy format not match, exit!")

        query_fp = q_BD == 'FP'
        query_tp = q_BD == 'TP'
        truth_fn = t_BD == 'FN'
        truth_tp = t_BD == 'TP'

        is_query_snp_fp = (q_BVT == 'SNP') and query_fp
        is_query_snp_tp = (q_BVT == 'SNP') and query_tp
        is_truth_snp_fn = (t_BVT == 'SNP') and truth_fn
        is_truth_snp_tp = (t_BVT == 'SNP') and truth_tp

        is_query_indel_fp = (q_BVT == 'INDEL') and query_fp
        is_query_indel_tp = (q_BVT == 'INDEL') and query_tp
        is_truth_indel_fn = (t_BVT == 'INDEL') and truth_fn
        is_truth_indel_tp = (t_BVT == 'INDEL') and truth_tp

        query_snp_fp = query_snp_fp + 1 if is_query_snp_fp else query_snp_fp
        query_snp_tp = query_snp_tp + 1 if is_query_snp_tp else query_snp_tp
        truth_snp_fn = truth_snp_fn + 1 if is_truth_snp_fn else truth_snp_fn
        truth_snp_tp = truth_snp_tp + 1 if is_truth_snp_tp else truth_snp_tp

        query_indel_fp = query_indel_fp + 1 if is_query_indel_fp else query_indel_fp
        query_indel_tp = query_indel_tp + 1 if is_query_indel_tp else query_indel_tp
        truth_indel_fn = truth_indel_fn + 1 if is_truth_indel_fn else truth_indel_fn
        truth_indel_tp = truth_indel_tp + 1 if is_truth_indel_tp else truth_indel_tp

        is_query_ins_fp = q_BI[0] == 'i' and is_query_indel_fp
        is_query_ins_tp = q_BI[0] == 'i' and is_query_indel_tp
        is_truth_ins_fn = t_BI[0] == 'i' and is_truth_indel_fn
        is_truth_ins_tp = t_BI[0] == 'i' and is_truth_indel_tp

        is_query_del_fp = q_BI[0] == 'd' and is_query_indel_fp
        is_query_del_tp = q_BI[0] == 'd' and is_query_indel_tp
        is_truth_del_fn = t_BI[0] == 'd' and is_truth_indel_fn
        is_truth_del_tp = t_BI[0] == 'd' and is_truth_indel_tp

        query_ins_fp = query_ins_fp + 1 if is_query_ins_fp else query_ins_fp
        query_ins_tp = query_ins_tp + 1 if is_query_ins_tp else query_ins_tp
        truth_ins_fn = truth_ins_fn + 1 if is_truth_ins_fn else truth_ins_fn
        truth_ins_tp = truth_ins_tp + 1 if is_truth_ins_tp else truth_ins_tp

        query_del_fp = query_del_fp + 1 if is_query_del_fp else query_del_fp
        query_del_tp = query_del_tp + 1 if is_query_del_tp else query_del_tp
        truth_del_fn = truth_del_fn + 1 if is_truth_del_fn else truth_del_fn
        truth_del_tp = truth_del_tp + 1 if is_truth_del_tp else truth_del_tp

    truth_all_tp = truth_snp_tp + truth_indel_tp
    truth_all_fn = truth_snp_fn + truth_indel_fn
    query_all_fp = query_snp_fp + query_indel_fp
    query_all_tp = query_snp_tp + query_indel_tp

    # p->precision, r->recall, f1->f1_score
    # a->overall, s->snp, id->indel, i->insertion, d->deletion
    ap, ar, af1 = metrics(query_fp=query_all_fp, query_tp=query_all_tp, truth_fn=truth_all_fn, truth_tp=truth_all_tp)
    sp, sr, sf1 = metrics(query_fp=query_snp_fp, query_tp=query_snp_tp, truth_fn=truth_snp_fn, truth_tp=truth_snp_tp)
    idp, idr, idf1 = metrics(query_fp=query_indel_fp, query_tp=query_indel_tp, truth_fn=truth_indel_fn, truth_tp=truth_indel_tp)
    ip, ir, if1 = metrics(query_fp=query_ins_fp, query_tp=query_ins_tp, truth_fn=truth_ins_fn, truth_tp=truth_ins_tp)
    dp, dr, df1 = metrics(query_fp=query_del_fp, query_tp=query_del_tp, truth_fn=truth_del_fn, truth_tp=truth_del_tp)

    print (''.join([item.ljust(20) for item in ["VariantType", 'TRUTH.FP', 'TRUTH.FN', 'TRUTH.TP','QUERY.TP', 'METRIC.Precision', 'METRIC.Recall', 'METRIC.F1_Score']]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["Overall", query_all_fp, truth_all_fn, truth_all_tp, query_all_tp, ap, ar, af1]]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["SNP", query_snp_fp, truth_snp_fn, truth_snp_tp, query_snp_tp, sp, sr, sf1]]),file=output_file)
    print (''.join([str(item).ljust(20) for item in ["INDEL", query_indel_fp, truth_indel_fn, truth_indel_tp, query_indel_tp, idp, idr, idf1]]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["INS", query_ins_fp, truth_ins_fn, truth_ins_tp, query_ins_tp, ip, ir, if1]]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["DEL", query_del_fp, truth_del_fn, truth_del_tp, query_del_tp, dp, dr, df1]]), file=output_file)
    print('\n', file=output_file)

    # print log_happy output
    pass_row = []
    snp_row = []
    indel_row = []
    if args.log_happy and os.path.exists(args.log_happy):
        log_happy = open(args.log_happy)
        for row in log_happy.readlines():
            if 'PASS' not in row:
                continue
            pass_row.append(row)

        for row in pass_row:

            if 'INDEL' in row:
                row = row.split()
                tp, fn, fp = row[3], row[4], row[6]
                precision, recall, f1 = row[11], row[10], row[13]
                indel_row = [fp, fn, tp, precision, recall, f1]
            if 'SNP' in row:
                row = row.split()
                tp, fn, fp = row[3], row[4], row[6]
                precision, recall, f1 = row[11], row[10], row[13]
                snp_row = [fp, fn, tp, precision, recall, f1]
        print('Double check with happy log:', file=output_file)
        print(' '.join(['%.6f' % item for item in [sp, sr, sf1, idp, idr, idf1]] + [str(item) for item in
                                                                                    [query_snp_fp, truth_snp_fn,
                                                                                     truth_snp_tp, query_indel_fp,
                                                                                     truth_indel_fn, truth_indel_tp]]),
              file=output_file)
        print(' '.join(snp_row[3:]) + ' ' + ' '.join(indel_row[3:]) + ' ' + ' '.join(snp_row[:3]) + ' ' + ' '.join(
            indel_row[:3]), file=output_file)
        print('\n', file=output_file)

    print(' '.join([str(item) for item in [ap, ar, af1, sp, sr, sf1, idp, idr, idf1, ip, ir, if1, dp, dr, df1]]),
          file=output_file)
    print(' '.join([str(item) for item in
                    [query_all_fp, truth_all_fn, truth_all_tp, query_snp_fp, truth_snp_fn, truth_snp_tp, query_indel_fp,
                     truth_indel_fn, truth_indel_tp, query_ins_tp, truth_ins_tp, truth_ins_tp, query_del_fp,
                     truth_del_fn, truth_del_tp]]), file=output_file)

    if output_fn:
        output_file.close()


def main():
    parser = ArgumentParser(description="Overall Metrics of hap.py output")

    parser.add_argument('--happy_vcf_fn', type=str, default=None,
                        help="Path to the happy vcf output file")

    parser.add_argument('--log_happy', type=str, default=None,
                        help="Path to the happy vcf output file")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Filename of the metrics output")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Cal(args)


if __name__ == "__main__":
    main()
