#!/usr/bin/env python

import argparse
import gzip
import os
import re
import sys
import time
from datetime import datetime
import pandas as pd
from Bio import SeqIO, bgzf


def main():
    start_time = time.time()
    now_time = datetime.today()
    parser = argparse.ArgumentParser("demultiplex PE fastq by paired reads classification of kraken2 output")
    parser.add_argument(
        '--r1',
        help='r1 fastq input'
    )
    parser.add_argument(
        '--r2',
        help='r2 fastq input'
    )
    parser.add_argument(
        '--kraken2_output',
        help='kraken2 output file'
    )
    parser.add_argument(
        '--taxonomy',
        default="genus",
        choices=["superkingdom", "phylum", "class", "order", "family", "genus"],
        help='highest taxonomy level, default: genus'
    )
    parser.add_argument(
        '--lineage',
        help='ncbi taxonomy lineage database'
    )
    parser.add_argument(
        '--outdir',
        help='output directory'
    )
    args = parser.parse_args()

    lineages = ["superkingdom", "phylum", "class", "order", "family", "genus"]
    lineages_ = ["unclassified"] + lineages

    if args.taxonomy not in lineages:
        print("%s not exists in ncbi taxonomy: [%s]" % (args.taxonomy, " ".join(lineages)))
        sys.exit(1)
    if not args.taxonomy:
        args.taxonomy = "genus"

    sub_lineages = lineages[0:(lineages.index(args.taxonomy)+1)]
    sub_lineages_ = lineages_[0:(lineages_.index(args.taxonomy)+1)]

    os.makedirs(args.outdir, exist_ok=True)
    for level in sub_lineages_:
        os.makedirs(os.path.join(args.outdir, level), exist_ok=True)

    log_h = open(os.path.join(args.outdir, "demultiplexer_kraken2.log"), 'w')

    log_h.write("now: %s\n" % now_time)
    taxonomy_lineages = pd.read_csv(args.lineage, sep='\t', low_memory=False)\
                          .set_index("tax_id", drop=False)
    
    log_h.write("step_1: parse taxonomy lineage has spent %d s\n" % (time.time() - start_time))

    taxonomys = {}
    demultiplexer = {}
    reads_partition = {}
    handle_r1 = {}
    handle_r2 = {}

    for level in sub_lineages_:
        taxonomys[level] = {}
        reads_partition[level] = {}
        handle_r1[level] = {}
        handle_r2[level] = {}

    taxonomys["unclassified"]["unclassified"] = 0

    with open(args.kraken2_output, 'r') as kh:
        for line in kh:
            read_id = line.split()[1]
            tax_id = int(line.split("(")[-1].split(')')[0].split()[-1])

            if tax_id == 0:
                taxonomys["unclassified"]["unclassified"] += 1
                demultiplexer[read_id] = ("unclassified", "unclassified")
            else:
                try:
                    lineage_dict = taxonomy_lineages.loc[tax_id, sub_lineages].to_dict()
                    have_taxonomy = False
                    for level in sub_lineages[::-1]:
                        taxonomy = lineage_dict[level]
                        if not pd.isnull(taxonomy):
                            taxonomy = "_".join(taxonomy.split())
                            if taxonomy not in taxonomys[level]:
                                taxonomys[level][taxonomy] = 1
                            else:
                                taxonomys[level][taxonomy] += 1
                            demultiplexer[read_id] = (level, taxonomy)
                            have_taxonomy = True
                            break
                    if not have_taxonomy:
                        taxonomys["unclassified"]["unclassified"] += 1
                        demultiplexer[read_id] = ("unclassified", "unclassified")
                        log_h.write("warning: taxid %d have no %s taxonomy level\n" % (tax_id, level))
                except KeyError:
                    taxonomys["unclassified"]["unclassified"] += 1
                    demultiplexer[read_id] = ("unclassified", "unclassified")

    log_h.write("step_2: parse kraken2_output has spent %d s, total %d taxnonmy level\n" % (time.time() - start_time, len(taxonomys)))

    total_reads_pair_1 = 0
    for level in sub_lineages_:
        out_dir = os.path.join(args.outdir, level)
        log_h.write("\t%s taxonomy level:\n" % level)
        for taxonomy in taxonomys[level]:
            log_h.write("\t\t%s reads num: %d\n" % (taxonomy, taxonomys[level][taxonomy]))
            total_reads_pair_1 += taxonomys[level][taxonomy]

            handle_r1[level][taxonomy] = bgzf.BgzfWriter(os.path.join(out_dir, taxonomy + ".1.fq.gz"), 'wb')
            handle_r2[level][taxonomy] = bgzf.BgzfWriter(os.path.join(out_dir, taxonomy + ".2.fq.gz"), 'wb')
        
            reads_partition[level][taxonomy] = {}
            reads_partition[level][taxonomy]["r1"] = []
            reads_partition[level][taxonomy]["r2"] = []

    if args.r1.endswith(".gz"):
        r1_handle = gzip.open(args.r1, 'rt')
        r2_handle = gzip.open(args.r2, 'rt')
    else:
        r1_handle = open(args.r1, 'r')
        r2_handle = open(args.r2, 'r')

    log_h.write("step_3: demultiplex paired reads\n")
    chunk_size = 100000
    count = 0
    total_reads_pair_2 = 0
    for read_1, read_2 in zip(SeqIO.parse(r1_handle, 'fastq'), SeqIO.parse(r2_handle, 'fastq')):
        total_reads_pair_2 += 1
        read_id = read_1.id.split("/")[0]
        if read_id in demultiplexer:
            level, taxonomy = demultiplexer[read_id]
            # print(level + " " + taxonomy)
            reads_partition[level][taxonomy]["r1"].append(read_1)
            reads_partition[level][taxonomy]["r2"].append(read_2)
            count += 1
        else:
            log_h.write("%s not in kraken2 output\n" % read_id)
        if count == chunk_size:
            for level in reads_partition:
                for taxonomy in reads_partition[level]:
                    if len(reads_partition[level][taxonomy]["r1"]) > 0:
                        SeqIO.write(reads_partition[level][taxonomy]["r1"], handle_r1[level][taxonomy], 'fastq')
                        SeqIO.write(reads_partition[level][taxonomy]["r2"], handle_r2[level][taxonomy], 'fastq')
                        reads_partition[level][taxonomy]["r1"] = []
                        reads_partition[level][taxonomy]["r2"] = []
            count = 0

    for level in reads_partition:
        for taxonomy in reads_partition[level]:
            if len(reads_partition[level][taxonomy]["r1"]) > 0:
                if taxonomy not in reads_partition[level]:
                    print("%s not in reads_partition 1" % taxonomy)
                    break
                if taxonomy not in reads_partition[level]:
                    print("%s not in reads_partition 2" % taxonomy)
                    break
                if taxonomy not in handle_r1[level]:
                    print("%s not in handle r1" % taxonomy)
                    break
                if taxonomy not in handle_r2[level]:
                    print("%s not in handle r2" % taxonomy)
                    break
                SeqIO.write(reads_partition[level][taxonomy]["r1"], handle_r1[level][taxonomy], 'fastq')
                SeqIO.write(reads_partition[level][taxonomy]["r2"], handle_r2[level][taxonomy], 'fastq')
                reads_partition[level][taxonomy]["r1"] = []
                reads_partition[level][taxonomy]["r2"] = []
            handle_r1[level][taxonomy].close()
            handle_r2[level][taxonomy].close()

    r1_handle.close()
    r2_handle.close()

    log_h.write("demultiplex %d paired reads has spent %s s\n" % (total_reads_pair_2, time.time() - start_time))
    log_h.write("total reads pair(kraken2): %d" % total_reads_pair_1)
    log_h.write("total reads pair(real): %d" % total_reads_pair_2)
    log_h.write("done\n")
    now_time = datetime.today()
    log_h.write("now: %s\n" % now_time)

    log_h.close()


if __name__ == '__main__':
    main()
