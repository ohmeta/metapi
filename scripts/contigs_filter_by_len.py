#!/usr/bin/env python
from contigs_to_gene import cut_fasta_by_len
import argparse

def main():
    parser = argparse.ArgumentParser(description="cut fasta by len")
    parser.add_argument('-fa', type=str, help='scaffolds or contigs file path')
    parser.add_argument('-sclen', type=int, help='scaffold or contigs length cutoff, default: 500', default=500)
    parser.add_argument('-outdir', type=str, help='output dir store gene prediction results')
    parser.add_argument('-prefix', type=str, help='prefix for file name')
    args = parser.parse_args()
    cut_fasta_by_len(args.fa, args.sclen, args.outdir, args.prefix, ".fa")

if __name__ == '__main__':
    main()