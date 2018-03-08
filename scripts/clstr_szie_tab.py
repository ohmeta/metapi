#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

def cdhit_clstr_parser(handle):
    """Generator function to iterate over cdhit clstr records (as string tuple)"""
    #Skip any text before the first record (e.g. blank lines, comments)
    #seq_id = ""
    #seq_len = ""
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in cdhit cluster file(fasta format) should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            # lines contain many cluster records
            lines.append(line.rstrip())
            line = handle.readline()

        # Remove trailing whitespace, and any internal spaces
        # (and any embedded \r which are possible in mangled files
        # when not opened in universal read lines mode)

        yield title, "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return  # StopIteration

    assert False, "Should not reach this line"

def clstr_size_tab(clstr_file, clstr_size_out):
    with open(clstr_file, 'r') as clstr_handle:
        for title, seq in cdhit_clstr_parser(clstr_handle):
            print(title + "\t" + seq)

def main():
    parser = argparse.ArgumentParser(description='parse cdhit cluster file and get cluste size distribution')
    parser.add_argument('--clstr', type=str, help='cluster file')
    parser.add_argument('--out', type=str, help='cluster size distribution')
    args = parser.parse_args()

    clstr_size_tab(args.clstr, args.out)

if '__name__' == '__main__':
    main()


