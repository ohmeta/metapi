#!/usr/bin/env python
import argparse
import gzip
from Bio import SeqIO, bgzf

def filter_pe_fasq_by_len(fq_1, fq_2, minlen, prefix):
    '''filter pe reads by min length'''
    fq_1_ = prefix + ".gt" + str(minlen) + ".1.fq.gz"
    fq_2_ = prefix + ".gt" + str(minlen) + ".2.fq.gz"
    with bgzf.BgzfWriter(fq_1_, 'wb') as out_1, bgzf.BgzfWriter(fq_2_, 'wb') as out_2:
        with gzip.open(fq_1, 'rt') as in_1, gzip.open(fq_2, 'rt') as in_2:
            for rec_a, rec_b in zip(SeqIO.parse(in_1, 'fastq'), SeqIO.parse(in_2, 'fastq')):
                if (len(rec_a.seq) > minlen) and (len(rec_b.seq) > minlen):
                    SeqIO.write(rec_a, out_1, 'fastq')
                    SeqIO.write(rec_b, out_2, 'fastq')

def main():
    '''main function'''
    parser = argparse.ArgumentParser(
        description='filter fastq file by reads length')
    parser.add_argument('-1', '--read1', help='paired-end fastq file one')
    parser.add_argument('-2', '--read2', help='paired-end fastq file two')
    parser.add_argument('-l', '--minlen', type=int, default=80,
                        help='remove reads if length < min-len')
    parser.add_argument('-p', '--prefix',
                        help='output prefix')
    args = parser.parse_args()

    filter_pe_fasq_by_len(args.read1, args.read2, args.minlen, args.prefix)

if __name__ == '__main__':
    main()
