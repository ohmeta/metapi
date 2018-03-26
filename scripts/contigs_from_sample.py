#!/usr/bin/env python
import sys
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser

def contigs_from_sample(contigs_len, sc_out):
    info = {}
    #count = 0
    with open(contigs_len, 'r') as handle:
        for line in handle:
            key = '_'.join(line.split("_")[:3])
            len = int(line.split("\t")[-1])
            if key not in info:
                info[key] = {}
                info[key]["num"] = 1
                info[key]["len"] = len
            else:
                info[key]["num"] += 1
                info[key]["len"] += len
            #count += 1
            #if count == 10000:
            #    break
    with open(sc_out, 'w') as out:
        out.write("sample_name\ttotal_contigs_num\ttotal_contigs_len\n")
        for key in info:
            out.write(key + "\t" + str(info[key]["num"]) + "\t" +
                      str(info[key]["len"]) + "\n")

def contigs_from_sample_list(contigs_list, sc_out):
    info = {}
    with open(contigs_list, 'r') as contigs_handle:
        for contigs_path in contigs_handle:
            key = os.path.basename(contigs_path.strip()).split(".")[0]
            if key not in info:
                info[key] = {}
                info[key]["num"] = 0
                info[key]["num_gt2kb"] = 0
                info[key]["len"] = 0
                info[key]["len_gt2kb"] = 0
            with open(contigs_path.strip(), 'r') as contigs_fa:
                for title, seq in SimpleFastaParser(contigs_fa):
                    info[key]["num"] += 1
                    info[key]["len"] += len(seq)
                    if len(seq) >= 2000:
                        info[key]["num_gt2kb"] += 1
                        info[key]["len_gt2kb"] += len(seq)
    with open(sc_out, 'w') as out:
        out.write("sample_name\ttotal_contigs_num\ttotal_contigs_num_gt2kb\ttotal_contigs_len\ttotal_contigs_len_gt2kb\n")
        for key in info:
            out.write("%s\t%d\t%d\t%d\t%d\n" % (key, info[key]["num"], info[key]["num_gt2kb"], info[key]["len"], info[key]["len_gt2kb"]))
            
def main():
    #contigs_from_sample(sys.argv[1], sys.argv[2])
    contigs_from_sample_list(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
