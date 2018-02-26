#!/usr/bin/env python
from Bio import SeqIO
import argparse


def merge_more2kb_contigs(contigs_list, merge_out):
    with open(merge_out, 'w') as out:
        with open(contigs_list, 'r') as contigs_handle:
            for contigs in contigs_handle:
                with open(contigs.strip(), 'r') as contigs_fa:
                    for record in  SeqIO.parse(contigs_fa,  'fasta'):
                        if len(record.seq) >= 2000:
                            SeqIO.write(record, out, 'fasta')
                        
def main():
    parser = argparse.ArgumentParser(description="merge more 2Kb contigs")
    parser.add_argument('--contigs_list', type=str, help='input file contain contigs path list')
    parser.add_argument('--merge_out', type=str, help='output file contain more2kb contigs')
    args = parser.parse_args()
    merge_more2kb_contigs(args.contigs_list, args.merge_out)


if __name__ == "__main__":
    main()