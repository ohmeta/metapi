#!/usr/bin/env python

from Bio.Alphabet import generic_dna
from Bio import SeqIO
import argparse
import re
import os
import time


def take_second(elem):
    return elem[1]


def parse_mgs(mgs_profile):
    bins = {}
    bins_num = []
    with open(mgs_profile, 'r') as ih:
        for line in ih:
            line_list = re.split(r"\s+|,", line.strip(",|\n"))
            bin_id = line_list[0]
            contigs_count = int(line_list[1])
            bins[bin_id] = []
            bins_num.append((bin_id, contigs_count))
            for contig_id in line_list[2:]:
                bins[bin_id].append(contig_id)
    bins_num.sort(key=take_second, reverse=True)
    return bins, bins_num


def extract(contigs_list, bins, bins_num, head, tail, outdir):
    files = []
    all_count = 0
    with open(contigs_list, 'r') as ih:
        for line in ih:
            files.append(line.strip())

    begin = time.time()
    records = SeqIO.index_db(":memory:", files, "fasta", generic_dna)
    end = time.time()
    print("index db: %.2f s" % (end - begin))

    begin = time.time()

    if head is not None:
        if head > len(bins_num):
            count = len(bins_num)
        else:
            count = head
        all_count += count
        for i in range(count):
            bin_id = bins_num[i][0]
            with open(os.path.join(outdir, bin_id + ".fa"), 'w') as oh:
                for contig_id in bins[bin_id]:
                    if contig_id in records:
                        SeqIO.write(records[contig_id], oh, 'fasta')
                    else:
                        print("%s has not find %s" % (bin_id, contig_id))

    if tail is not None:
        if tail > len(bins_num):
            count = len(bins_num)
        else:
            count = tail
        all_count += count
        for i in range(count):
            bin_id = bins_num[-(i+1)][0]
            with open(os.path.join(outdir, bin_id + ".fa"), 'w') as oh:
                for contig_id in bins[bin_id]:
                    if contig_id in records:
                        SeqIO.write(records[contig_id], oh, 'fasta')
                    else:
                        print("%s has not find %s" % (bin_id, contig_id))

    if (head is None) and (tail is None):
        for bin_id in bins:
            all_count += 1
            with open(os.path.join(outdir, bin_id + ".fa"), 'w') as oh:
                for contig_id in bins[bin_id]:
                    if contig_id in records:
                        SeqIO.write(records[contig_id], oh, 'fasta')
                    else:
                        print("%s has not find %s" % (bin_id, contig_id))

    records.close()
    end = time.time()

    print("extract %d bins: %.2f s" % (all_count, end - begin))


def main():
    parser = argparse.ArgumentParser("get bins fasta from mgs contigs/scafoolds profile")
    parser.add_argument('-p', '--profile', type=str, help='mgs contigs/scaffolds profile')
    parser.add_argument('-l', '--contigs_list', type=str, help='assembly contigs/scaffolds fasta path list')
    parser.add_argument('-o', '--outdir', type=str, help='bins output dir')
    parser.add_argument('--head', type=int, default=None, help='head number bins')
    parser.add_argument('--tail', type=int, default=None, help='tail number bins')

    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    (bins, bins_num) = parse_mgs(args.profile)

    if (args.head is not None) and (args.tail is not None):
        assert args.head + args.tail <= len(bins_num), "too many head or too many tail"

    extract(args.contigs_list, bins, bins_num, args.head, args.tail, args.outdir)


if __name__ == '__main__':
    main()
