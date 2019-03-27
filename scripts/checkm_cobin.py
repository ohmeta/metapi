#!/usr/bin/env python

from Bio import SeqIO
import argparse
import re

def parse_mgs(mgs_profile):
    cag = {}
    with open(mgs_profile, 'r') as ih:
        for line in ih:
            line_list = re.split(r"\s+|,", line)
            cag_id = line_list[0]
            # seq_count = line_list[1]
            cag[cag_id] = {}
            for contig_id in line_list[2:]:
                sample_id = contig_id.split('_')[0]
                if sample_id not in cag[cag_id]:
                    cag[cag_id][sample_id] = []
                if "cov" in contig_id:
                    spades_contig_id = "node_" + "_".join(contig_id.split("_")[1:])
                    cag[cag_id][sample_id].append(spades_contig_id)
                else:
                    megahit_contig_id = "K199" + "_".join(contig_id.split("_")[1:])
                    cag[cag_id][sample_id].append(megahit_contig_id)
    return cag


def contigs_index():
    pass


def main():
    parser = argparse.ArgumentParser("get bins fasta from mgs contigs/scafoolds profile")
    parser.add_argument('-p', '--profile', type=str, help='mgs contigs/scaffolds profile')
    parser.add_argument('-l', '--contigs_list', type=str, help='assembly contigs/scaffolds fasta path list')
    parser.add_argument('-o', '--outdir', type=str, help='bins output dir')

    args = parser.parse_args()



if __name__ == '__main__':
    main()