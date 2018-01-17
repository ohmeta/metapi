#!/usr/bin/env python
## Metagenomics Institute of BGI Research
## zhujie@genomics.cn
## 2017-12-05
## GPL-V3

import subprocess
import argparse
import os
import csv
import gzip
from Bio import SeqIO

def gen_statout_matrix(fqfile_1, fqfile_2, fqfile_s, statout_f, minl, maxl):
    allfq_info_dict = {}
    allfq_info_dict['header'] = []
    allfq_info_dict['body'] = []
    gc_qual_l80_key = ['GC', 'Q10', 'Q20', 'Q30',
                       'total_base_num', 'total_reads_num',
                       'L80_base_num', 'L80_reads_num',
                       'L80_base_rate', 'L80_reads_rate']
    for i in gc_qual_l80_key:
        for j in ['', '_1', '_2', '_s']:
            allfq_info_dict['header'].append(i + j)
    for i in ['bp', 'bp_1', 'bp_2', 'bp_s']:
        for j in range(minl, maxl + 1):
            allfq_info_dict['header'].append(str(j) + i)

    statout_key = ['sample_name',  'host_rate',  'all_align_reads_num',
                   'Total_pair',   'Aligned_pair',
                   'Total_single', 'Aligned_single',
                   'read1_base',   'read2_base', 'single_base',
                   'max_length',   'min_length', 'avg_length']

    allfq_info_dict["header"] = statout_key + allfq_info_dict["header"]

    sample_name = os.path.basename(statout_f).split('.')[0]
    stat_out_dict = parse_statout(statout_f, sample_name)

    len_gc_qual = {}
    len_gc_qual_1 = get_len_gc_qual(fqfile_1)
    len_gc_qual_2 = get_len_gc_qual(fqfile_2)
    len_gc_qual_s = get_len_gc_qual(fqfile_s)

    for i in range(minl, maxl + 1):
        key_1 = str(i) + get_suffix(fqfile_1)
        key_2 = str(i) + get_suffix(fqfile_2)
        key_s = str(i) + get_suffix(fqfile_s)
        key   = str(i) + get_suffix(statout_f)
        if key_1 not in len_gc_qual_1:
            len_gc_qual_1[key_1] = 0
        if key_2 not in len_gc_qual_2:
            len_gc_qual_2[key_2] = 0
        if key_s not in len_gc_qual_s:
            len_gc_qual_s[key_s] = 0
        len_gc_qual[key] = len_gc_qual_1[key_1] + \
                           len_gc_qual_2[key_2] + \
                           len_gc_qual_s[key_s]

    # len_gc_qual['total_base_num]  == stat_out_dict['reads1_base'] + \
    #                                  stat_out_dict['reads2_base'] + \
    #                                  stat_out_dict['single_base']
    # len_gc_qual['total_reads_num] == stat_out_dict['Aligned_pair'] * 2 + \
    #                                  stat_out_dict['Aligned_single']
    for key in gc_qual_l80_key[4:8]:
        len_gc_qual[key] = len_gc_qual_1[key + '_1'] + \
                           len_gc_qual_2[key + '_2'] + \
                           len_gc_qual_s[key + '_s']
    len_gc_qual["L80_base_rate"]  = len_gc_qual['L80_base_num']  / len_gc_qual['total_base_num']
    len_gc_qual["L80_reads_rate"] = len_gc_qual['L80_reads_num'] / len_gc_qual['total_reads_num']

    for key in gc_qual_l80_key[0:4]:
        len_gc_qual[key] = (len_gc_qual_1[key + '_1'] * len_gc_qual_1['total_base_num_1'] +
                            len_gc_qual_2[key + '_2'] * len_gc_qual_2['total_base_num_2'] +
                            len_gc_qual_s[key + '_s'] * len_gc_qual_s['total_base_num_s']) / \
                            len_gc_qual['total_base_num']
    ## merge many dicts to a dict
    fq_info_dict = { **stat_out_dict, **len_gc_qual, **len_gc_qual_1, **len_gc_qual_2, **len_gc_qual_s }
    allfq_info_dict['body'].append(fq_info_dict)

    return allfq_info_dict

def parse_statout(statout_f, sample_name):
    with open(statout_f, 'r') as h:
        stat_out_dict = { k : v
                          for k,v in zip(h.readline().strip().split(),
                                         h.readline().strip().split()) }
    stat_out_dict["sample_name"] = sample_name
    all_align_reads_num = (int(stat_out_dict["Total_pair"]) - int(stat_out_dict["Aligned_pair"])) * 2 + \
                          int(stat_out_dict["Total_single"]) - int(stat_out_dict["Aligned_single"])
    host_rate = all_align_reads_num / \
                (int(stat_out_dict["Total_pair"]) * 2 + int(stat_out_dict["Total_single"]))
    stat_out_dict["all_align_reads_num"] = all_align_reads_num
    stat_out_dict["host_rate"] = host_rate

    return stat_out_dict

def lentab(fqfile):
    ## https://github.com/alastair-droop/fqtools
    len_tab_dict = {}
    suffix       = get_suffix(fqfile)
    L80_key      = "L80_rate" + suffix.lstrip('bp')
    L80_count    = 0
    total_count  = 0
    cmd          = "fqtools lengthtab " + fqfile
    len_tab_seq = subprocess.getoutput(cmd)
    for d in len_tab_seq.strip().split("\n"):
        one_row = d.strip().split('\t')
        key = str(one_row[0]) + suffix
        len_tab_dict[key] = int(one_row[1])
        total_count += int(one_row[1])
        if int(one_row[0]) >= 80:
            L80_count += int(one_row[1])
    len_tab_dict[L80_key] = L80_count / total_count

    return len_tab_dict

def get_len_gc_qual(fqfile):
    len_gc_qual_dict = {}
    suffix_          = get_suffix(fqfile)
    suffix           = suffix_.lstrip("bp")
    gc_num           = 0
    q10_num          = 0
    q20_num          = 0
    q30_num          = 0
    total_base_num   = 0
    total_reads_num  = 0
    l80_base_num     = 0
    l80_reads_num    = 0
    with gzip.open(fqfile, 'rt') as fq:
        for rec in SeqIO.parse(fq, 'fastq'):
            seq_len          = len(rec.seq)
            total_reads_num += 1
            total_base_num  += seq_len
            len_key          = str(seq_len) + suffix_
            if seq_len >= 80:
                l80_reads_num += 1
                l80_base_num  += seq_len
            if len_key in len_gc_qual_dict:
                len_gc_qual_dict[len_key] += 1
            else:
                len_gc_qual_dict[len_key] = 1
            gc_num += rec.seq.count('G') + rec.seq.count('C')
            for qual in rec.letter_annotations['phred_quality']:
                if qual >= 30:
                    q30_num += 1
                    q20_num += 1
                    q10_num += 1
                elif qual >= 20:
                    q20_num += 1
                    q10_num += 1
                elif qual >= 10:
                    q10_num += 1
    len_gc_qual_dict['total_base_num'  + suffix]  = total_base_num
    len_gc_qual_dict['total_reads_num' + suffix]  = total_reads_num
    len_gc_qual_dict['L80_base_num'    + suffix]  = l80_base_num
    len_gc_qual_dict['L80_reads_num'   + suffix]  = l80_reads_num
    len_gc_qual_dict['L80_base_rate'   + suffix]  = l80_base_num  / total_base_num
    len_gc_qual_dict['L80_reads_rate'  + suffix]  = l80_reads_num / total_reads_num
    len_gc_qual_dict['GC'              + suffix]  = gc_num        / total_base_num
    len_gc_qual_dict['Q10'             + suffix]  = q10_num       / total_base_num
    len_gc_qual_dict['Q20'             + suffix]  = q20_num       / total_base_num
    len_gc_qual_dict['Q30'             + suffix]  = q30_num       / total_base_num

    return len_gc_qual_dict

def get_suffix(fqfile):
    if fqfile.endswith("1.fq.gz"):
        suffix = "bp_1"
    elif fqfile.endswith("2.fq.gz"):
        suffix = "bp_2"
    elif fqfile.endswith("single.fq.gz"):
        suffix = "bp_s"
    else:
        suffix = "bp"

    return suffix

def write_csv(outdir, statout_f, header, body):
    sample_name = os.path.basename(statout_f).split('.')[0]
    outfile = os.path.join(outdir, sample_name + ".tsv")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    with open(outfile, 'w') as tsvfile:
        tsvwriter = csv.DictWriter(tsvfile, header, delimiter = '\t')
        tsvwriter.writeheader()
        tsvwriter.writerows(body)

def main():
    parser = argparse.ArgumentParser(description='generate rmhost fastq sequence length distribution and stat_out matrix')
    parser.add_argument('-1', dest='fq1',            required=True, help='xx.1.fq.gz')
    parser.add_argument('-2', dest='fq2',            required=True, help='xx.2.fq.gz')
    parser.add_argument('-s', dest='fqs',            required=True, help='xx.single.fq.gz')
    parser.add_argument('-r', dest='statout',        required=True, help='xx.stat_out')
    parser.add_argument('-n', dest='minl', type=int, required=True, help='reads min length')
    parser.add_argument('-m', dest='maxl', type=int, required=True, help='reads max length')
    parser.add_argument('-o', dest='outdir',         required=True, help='a dir store matrix output csv formated')
    args = parser.parse_args()
    allfq_info_dict = gen_statout_matrix(args.fq1, args.fq2, args.fqs, args.statout, args.minl, args.maxl)
    write_csv(args.outdir, args.statout, allfq_info_dict["header"], allfq_info_dict["body"])

if __name__ == "__main__":
    main()
