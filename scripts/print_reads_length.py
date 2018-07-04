#!/usr/bin/env python
"get each reads length form fasta/fastq file"
import argparse
import gzip
from Bio import SeqIO

def print_len(infile, seqtype):
    '''print_len function'''
    if infile.endswith(".gz"):
        handle = gzip.open(infile, 'rt')
    else:
        handle = open(infile, 'rt')
    for reads in SeqIO.parse(handle, seqtype):
        print(reads.id, "\t", len(reads))
    handle.close()


def main():
    '''main function'''
    parser = argparse.ArgumentParser(description='print each reads id and length info form fasta/fastq file')
    parser.add_argument('--infile', action='store', dest='infile', help='input fasta/fastq file')
    parser.add_argument('--seqtype', action='store', dest='seqtype', help='input file seq type, fasta or fastq')
    args = parser.parse_args()
    print_len(args.infile, args.seqtype)

if __name__ == "__main__":
    main()
