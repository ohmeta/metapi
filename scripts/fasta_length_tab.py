#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
"""get each sequence length from a fasta file and pring it to a file, then plot"""

def gen_fa_len_tab(fa_file, len_out):
    with open(len_out, 'w') as out_handle:
        with open(fa_file, 'r') as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                #out_handle.write(title + "\t" + str(len(seq)))
                # just print id and seq length
                out_handle.write(title.split(' ')[0] + "\t" + str(len(seq)) + "\n")

# megahit contigs header contains contigs length info
def gen_fa_len_tab_megahit(fa_file, len_out):
    with open(len_out, 'w') as out_handle:
        with open(fa_file, 'r') as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                # maybe wrong
                len = title.split(' ')[-1].split('=')[-1]
                out_handle.write(title + "\t" + len + "\n")

def main():
    parser = argparse.ArgumentParser(description='get fasta length info')
    parser.add_argument('--fasta', type=str, help='fasta file')
    parser.add_argument('--out', type=str, help='fasta length output file')
    args = parser.parse_args()

    gen_fa_len_tab(args.fasta, args.out)
    
    # fasta input must contigs file which was assemblyed by megahit
    # gen_fa_len_tab(args.fasta, args.out)

if __name__ == '__main__':
    main()