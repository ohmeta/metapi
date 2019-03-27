#!/usr/bin/env python

from Bio.Alphabet import generic_dna
from Bio import SeqIO
import argparse
import re
import os


def parse_mgs(mgs_profile):
    bins = {}
    with open(mgs_profile, 'r') as ih:
        for line in ih:
            line_list = re.split(r"\s+|,", line)
            bin_id = line_list[0]
            bins[bin_id] = []
            for contig_id in line_list[2:]:
                bins[bin_id].append(contig_id)
    return bins


def extract(contigs_list, bins, outdir):
    files = []
    for i in contigs_list:
        files.append(i)
    records = SeqIO.index_db(":memory:", files, "fasta", generic_dna)
    for bin_id in bins:
        with open(os.path.join(outdir, "bin_id.fa"), 'w') as oh:
            for contig_id in bins[bin_id]:
                SeqIO.write(records[contig_id], oh, 'fasta')


def main():
    parser = argparse.ArgumentParser("get bins fasta from mgs contigs/scafoolds profile")
    parser.add_argument('-p', '--profile', type=str, help='mgs contigs/scaffolds profile')
    parser.add_argument('-l', '--contigs_list', type=str, help='assembly contigs/scaffolds fasta path list')
    parser.add_argument('-o', '--outdir', type=str, help='bins output dir')

    args = parser.parse_args()
    bins = parse_mgs(args.profile)
    extract(args.contigs_list, bins, args.outdir)


if __name__ == '__main__':
    main()
