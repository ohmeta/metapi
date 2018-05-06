#!/usr/bin/env python
import argparse

def aggregate(cov):
    '''
    bedtools genomecov -ibam sample.mapped.sorted.bam -g contigs_c10K.len > sample_cov.txt
    produce a histogram of coverage of the exons throughout the genome
    
    output format explain:
    1. chromosome(or entire genome)
    2. depth of coverage from features in input file
    3. number of bases on chromosome(or genome) with depth equal to column 2
    4. size of chromosome(or entire genome) in base pairs
    5. fraction of baes on chromosome(or entire genome) with depth equal to column 2
    so column_5 = column_3 / column_4
    all sum(column_3{column_1}) = column_4{column_1}
    all sum(column_5{column_1}) = 1

    k119_2  1       30      399     0.075188
    k119_2  2       27      399     0.0676692
    k119_2  3       151     399     0.378446
    k119_2  4       79      399     0.197995
    k119_2  5       54      399     0.135338
    k119_2  6       39      399     0.0977444
    k119_2  7       19      399     0.047619
    k119_3  0       387     473     0.818182
    k119_3  1       86      473     0.181818
    k119_4  4       1       340     0.00294118
    '''
    with open(cov, 'r') as in_handle:
        cov_num = {}
        chr_len = {}
        chr_list = []
        for line in in_handle:
            chr, depth, num, len, frac = line.strip().split('\t')
            if chr not in chr_len:
                chr_len[chr] = int(len)
                cov_num[chr] = int(depth) * int(num)
                chr_list.append(chr)
            else:
                cov_num[chr] += int(depth) * int(num)
        for chr_name in chr_list:
            print("%s,%f" % (chr_name, cov_num[chr_name] / chr_len[chr_name]))

def main():
    parser = argparse.ArgumentParser(description='aggregate the output of bedtools')
    parser.add_argument('-cov', type=str, help='input coverage file')
    args = parser.parse_args()

    aggregate(args.cov)

if __name__ == '__main__':
    main()

# awk -F'\t' '{l[$1]=l[$1]+($2*$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}'