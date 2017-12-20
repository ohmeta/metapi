#!/usr/bin/env python
## Metagenomics Institute of BGI Research
## zhujie@genomics.cn
## 2017-11-29
## GPL-V3

import os
import argparse

##TODO
## clean and SE

def parse_pe_clean_statout(handle, min_l, max_l):
    header_list = [str(i) for i in range(min_l, max_l + 1)]
    value_dict = {}
    for key in header_list:
        value_dict[key] = 0

    ## total info
    for key,value in zip(handle.readline().strip().split(), handle.readline().strip().split()):
        header_list.append(key)
        value_dict[key] = value

    ## reads_1 info
    tag = True
    for key,value in zip(handle.readline().strip().split(), handle.readline().strip().split()):
        if tag:
            tag = False
        else:
            key = key + "_1"
        header_list.append(key)
        value_dict[key] = value

    ## reads_2 info
    tag = True
    for key,value in zip(handle.readline().strip().split(), handle.readline().strip().split()):
        if tag:
            tag = False
        else:
            key = key + "_2"
        header_list.append(key)
        value_dict[key] = value

    ## reads_single info
    tag = True
    for key,value in zip(handle.readline().strip().split(), handle.readline().strip().split()):
        if tag:
            tag = False
        else:
            key = key + "_single"
        header_list.append(key)
        value_dict[key] = value

    ## length info
    next(handle)
    total_filter_base = 0
    total_filter_reads = 0
    total_filter_reads_len_gt80 = 0
    L80 = 0
    for line in handle:
        line_list = line.strip().split()
        reads_len = line_list[0]
        reads_num = line_list[1]
        total_filter_reads += int(reads_num)
        total_filter_base  += int(reads_len) * int(reads_num)
        if (int(reads_len) >= 80):
            total_filter_reads_len_gt80 += int(reads_num)
        value_dict[str(reads_len)] = str(reads_num)

    L80 = total_filter_reads_len_gt80 / total_filter_reads
    header_list.append("total_filter_base")
    value_dict["total_filter_base"] = str(total_filter_base)
    header_list.append("total_filter_reads")
    value_dict["total_filter_reads"] = str(total_filter_reads)
    header_list.append("L80")
    value_dict["L80"] = str(L80)

    return (value_dict, header_list)

def gen_len_matrix(dirname, min_l, max_l):
    no_header = True
    for fl in os.listdir(dirname):
        if fl.endswith("stat_out"):
            sample_name = fl.split("/")[-1].split(".")[0]
            statout = os.path.join(dirname, fl)
            with open(statout, 'r') as h:
                tuple_ = parse_pe_clean_statout(h, min_l, max_l)
                if no_header:
                    header = "sample_name\t" + "\t".join(tuple_[1])
                    print(header)
                    body = sample_name + "\t" + "\t".join([tuple_[0][key] for key in tuple_[1]])
                    print(body)
                    no_header = False
                else:
                    body = sample_name + "\t" + "\t".join([tuple_[0][key] for key in tuple_[1]])
                    print(body)

def main():
    parser = argparse.ArgumentParser("convert many clean statout to a matrix\n \
                                      e.g: python clean_statout_to_matrix.py ../data/clean_statout -m 100 -n 30 > ../data/length_clean_statout.tsv\n")
    parser.add_argument("-d", "--statout_dir", help="a directory contain many samples clean statout file")
    parser.add_argument("-m", "--max_len", type=int, help="max reads length")
    parser.add_argument("-n", "--min_len", type=int, help="min reads length")
    args = parser.parse_args()
    gen_len_matrix(args.statout_dir, args.min_len, args.max_len)

if __name__ == "__main__":
    main()
