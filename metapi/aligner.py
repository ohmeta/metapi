#!/usr/bin/env python
import argparse
import csv
import os
import re
from decimal import *


"""
How to assess the quality of metagenomcis assembly
https://www.biostars.org/p/128629/#128639

Brian Bushnell saied:
calculate the percentage of reads that map back to the assembly
if only 50% of your reads map to the assembly, it is not very complete
but if 95% of your reads map to the assembly, then even if it is
somewhat fragmented, that's probably very good

It might also help to look at the percentage properly paired reads to
detect any chimeras, something that seems especially relevant in a
metagenome assembly

the most useful tools is quast, quality assessment tool for
genome ascalculate map rate from bamsemblies

http://genomebio.org/alignment-stats-bwa/
getting alignment stats out of bwa

bwa mem -t 6 ref read.1.fq read.2.fq \
| samtools view -@6 -Sbh - \
| tee >(samtools flagstat - > stats.out) > aln.bam

http://www.pnas.org/content/pnas/113/42/11901.full.pdf
Deep sequencing of 10,000 human genomes(Amalio Telenti and J.Craig Venter)

metabat_coverage
concoct_coverage
checkm_coverage
"""


def flagstats_summary(flagstats, out_file, method):
    """
    get alignment rate from sorted bam file
    samtools -flagstat --threads 8 sample.sort.bam
    """
    headers = [
        "sample_id",
        "total_num",
        "read_1_num",
        "read_2_num",
        "mapping_type",
        "mapped_num",
        "mapped_rate",
        "paired_num",
        "paired_rate",
        "singletons_num",
        "singletons_rate",
        "mate_mapped_num",
        "mate_mapped_num_mapQge5",
    ]
    mapping_info = []
    getcontext().prec = 8

    # with open(flagstat_list, 'r') as list_handle:
    if method == 1:
        list_handle = open(flagstats, "r")
    if method == 2:
        list_handle = flagstats

    for flagstat_file in list_handle:
        if os.path.exists(flagstat_file.strip()):
            info = {}
            info["sample_id"] = os.path.basename(flagstat_file.strip()).split(".")[0]
            stat_list = open(flagstat_file.strip(), "r").readlines()
            info["total_num"] = stat_list[0].split(" ")[0]
            info["read_1_num"] = stat_list[6].split(" ")[0]
            info["read_2_num"] = stat_list[7].split(" ")[0]

            mapped = re.split(r"\(|\s+", stat_list[4])
            info["mapped_num"] = mapped[0]
            info["mapped_rate"] = Decimal(mapped[5].rstrip("%")) / Decimal(100)

            paired = re.split(r"\(|\s+", stat_list[8])
            info["paired_num"] = paired[0]
            paired_rate = paired[6].rstrip("%")
            if paired_rate != "N/A":
                info["paired_rate"] = Decimal(paired_rate) / Decimal(100)
                info["mapping_type"] = "paired-end"
            else:
                info["paired_rate"] = paired_rate
                info["mapping_type"] = "single-end"

            singletons = re.split(r"\(|\s+", stat_list[-3])
            info["singletons_num"] = singletons[0]
            singletons_rate = singletons[5].rstrip("%")
            if singletons_rate != "N/A":
                info["singletons_rate"] = Decimal(singletons_rate) / Decimal(100)
            else:
                info["singletons_rate"] = singletons_rate

            info["mate_mapped_num"] = re.split(r"\(|\s+", stat_list[-2])[0]
            info["mate_mapped_num_mapQge5"] = re.split(r"\(|\s+", stat_list[-1])[0]
            mapping_info.append(info)

    with open(out_file, "w") as out_handle:
        f_tsv = csv.DictWriter(out_handle, headers, delimiter="\t")
        f_tsv.writeheader()
        f_tsv.writerows(mapping_info)


def main():
    """main funciton"""
    parser = argparse.ArgumentParser(
        description="compute alignment rate based bam file"
    )
    parser.add_argument(
        "-statlist", default=None, type=str, help="sorted bam file list"
    )
    parser.add_argument("-statfiles", default=None, nargs="*", help="sorted bam file")
    parser.add_argument("-outfile", type=str, help="output alignment rate file")
    args = parser.parse_args()
    if args.statlist:
        method = 1
        flagstats_summary(args.statlist, args.outfile, method)
    if args.statfiles:
        method = 2
        flagstats_summary(args.statfiles, args.outfile, method)


if __name__ == "__main__":
    main()
