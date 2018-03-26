#!/usr/bin/env python
#from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import re

pattern = re.compile(r'\d+\t(\d+)[a-z]{2}, >(.+)\.\.\. \*')
#pattern = re.compile(r'\d+\t(\d+)[a-z]{2},\s>(.+)\.\.\.\s\*')
#pattern = re.compile(r'\d+\t(\d+)nt, >(.+)\.\.\. \*')
#pattern = re.compile(r'\d+\t(\d+)nt,\s>(.+)\.\.\.\s\*')

# this parser base code comes from Bio.SeqIO.FastaIO.SimpleFastaParser :)
def cdhit_clstr_parser(handle):
    """Generator function to iterate over cdhit clstr records (as string tuple)
    
    >Cluster 0
    0       1131322nt, >k119_12676... *
    1       84315nt, >k119_210239... at -/99.66%
    2       73592nt, >k119_187067... at +/99.86%
    3       70665nt, >k119_160147... at -/99.32%
    4       66352nt, >k119_217379... at +/99.89%
    5       63337nt, >k119_125106... at +/99.28%
    6       63232nt, >k119_150147... at -/99.80%
    7       59840nt, >k119_197728... at +/99.04%
    8       59306nt, >k119_59391... at -/99.00%
    >Cluster 5343379
    0       2000nt, >k119_192744... *
    >Cluster 5343380
    0       2000nt, >k119_222307... *
    >Cluster 5343381
    0       2000nt, >k119_232332... *
    >Cluster 5343382
    0       2000nt, >k119_241124... *
    >Cluster 5343383
    0       2000nt, >k119_253638... *
    
    """
    #Skip any text before the first record (e.g. blank lines, comments)
    seq_id = ""
    seq_len = 0
    clstr_size = 0
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        clstr_size = 0
        if line[0] != ">":
            raise ValueError(
                "Records in cdhit cluster file(fasta format) should start with '>' character")
        clstr_name = line[1:].rstrip()
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            # lines contain many cluster records
            #lines.append(line.rstrip())
            clstr_size += 1
            matches = re.search(pattern, line)
            if matches:
                seq_len = matches.group(1)
                seq_id = matches.group(2)

            line = handle.readline()
        yield clstr_name, seq_id, seq_len, clstr_size

        if not line:
            return  # StopIteration

    assert False, "Should not reach this line"

def clstr_size_tab(clstr_file, clstr_size_out):
    with open(clstr_size_out, 'w') as out_handle:
        out_handle.write("cluster_name\tcluster_size\tsequence_id\tsequence_length\n")
        with open(clstr_file, 'r') as clstr_handle:
            for clstr_name, seq_id, seq_len, clstr_size in cdhit_clstr_parser(clstr_handle):
                clstr_name = "cluster_" + clstr_name.split(' ')[1]
                out_handle.write(clstr_name + "\t" + str(clstr_size) + "\t" + seq_id + "\t" + str(seq_len) + "\n")

def main():
    parser = argparse.ArgumentParser(description='parse cdhit cluster file and get cluste size distribution')
    parser.add_argument('--clstr', type=str, help='cluster file')
    parser.add_argument('--out', type=str, help='cluster size distribution')
    args = parser.parse_args()

    clstr_size_tab(args.clstr, args.out)

if __name__ == '__main__':
    main()


