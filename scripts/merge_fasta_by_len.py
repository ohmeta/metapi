#!/usr/bin/env python
import argparse
import gzip
import os

from Bio import SeqIO


def merge_fa_by_len(falist, minlen, maxlen, outfa):
    with open(falist, 'r') as falist_h, open(outfa, 'w') as out_h:
        for fa_file in falist_h:
            fa_file = fa_file.rstrip()
            if fa_file.endswith(".gz"):
                fa_h = gzip.open(fa_file, 'rt')
            else:
                fa_h = open(fa_file, 'r')
            for record in SeqIO.parse(fa_h, 'fasta'):
                if (len(record.seq) >= minlen) and (len(record.seq) <= maxlen):
                    SeqIO.write(record, out_h, 'fasta')
            fa_h.close()


def main():
    parser = argparse.ArgumentParser(description="merge many fasta file to a fasta file by length cutoff")
    parser.add_argument('--falist', type=str, help='input file contain fasta path list')
    parser.add_argument('--minlen', type=int, help='sequences min length cutoff', default=1)
    parser.add_argument('--maxlen', type=int, help='sequences max length cutoff', default=10000000000)
    parser.add_argument('--outfa', type=str, help='output fasta contain sequences which length between [minlen, maxlen]')
    args = parser.parse_args()

    merge_fa_by_len(args.falist, args.minlen, args.maxlen, args.outfa)


if __name__ == "__main__":
    main()
