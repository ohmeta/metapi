#!/usr/bin/env python
import argparse
import os
import csv

def gen_size_tsv(fqlist, ctglist, tsvout):
    '''gen data size tsv out'''
    fq_size = {}
    ctg_size = {}
    file_size = {}
    file_size["header"] = ["fq_1", "fq_2", "fq_s", "contig", "sample_name"]
    file_size["body"] = []

    with open(fqlist, 'r') as fq_handle, open(ctglist, 'r') as ctg_handle:
        for (fq_line, ctg_line) in zip(fq_handle, ctg_handle):
            (reads_a, reads_b, reads_s) = fq_line.strip().split()
            fq_name = os.path.basename(reads_a).split('.')[0]
            ctg_name = os.path.basename(ctg_line).split('.')[0]
            if fq_name not in fq_size:
                fq_size[fq_name] = {}
            fq_size[fq_name]["fq_1"] = os.path.getsize(reads_a)
            fq_size[fq_name]["fq_2"] = os.path.getsize(reads_b)
            fq_size[fq_name]["fq_s"] = os.path.getsize(reads_s)
            if ctg_name not in ctg_size:
                ctg_size[ctg_name] = {}
            ctg_size[ctg_name] = os.path.getsize(ctg_line.strip())

    assert sorted(fq_size.keys()) == sorted(ctg_size.keys())

    for key in ctg_size:
        file_size_ = {}
        file_size_["sample_name"] = key
        file_size_["fq_1"] = fq_size[key]["fq_1"]
        file_size_["fq_2"] = fq_size[key]["fq_2"]
        file_size_["fq_s"] = fq_size[key]["fq_s"]
        file_size_["contig"] = ctg_size[key]
        file_size["body"].append(file_size_)

    with open(tsvout, 'w') as out_handle:
        f_tsv = csv.DictWriter(out_handle, file_size["header"], delimiter='\t')
        f_tsv.writeheader()
        f_tsv.writerows(file_size["body"])


def main():
    '''main function'''
    parser = argparse.ArgumentParser(
        description='''research relationships between fastq size and contigs size:
Usage: python fastq_contig_size_relationship.py --fqlist ./212S_rmhost_fqgz.pathlist.paired --ctglist ./212S_assembly_contigs.pathlist --tsvout fq_contigs_size.ts
''')
    parser.add_argument('--fqlist', type=str,
                        help='rmhost fastq file path list')
    parser.add_argument('--ctglist', type=str,
                        help='contigs file path list')
    parser.add_argument('--tsvout', type=str,
                        help='tsv out put')
    args = parser.parse_args()
    gen_size_tsv(args.fqlist, args.ctglist, args.tsvout)


if __name__ == '__main__':
    main()
