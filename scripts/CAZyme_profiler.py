#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
from xopen import xopen
import math
import concurrent.futures
import os

# example
# seqtk mergepe s1.1.fastq.gz s1.2.fastq.gz | \
#    diamond blastx -d CAZyDB.09242021.dmnd \
#    -o blastx_output/s1.m6 -f 6 -p 24 -e 1e-05 > s1.log 2>&1

# BLAST tab format
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# qseqid sseqid pident length evalue


def parse_fai(fai):
    df = pd.read_csv(fai, sep="\t", header=None,
    names=["CAZy_gene_name", "CAZy_gene_len", "CAZy_gene_start", "CAZy_gene_end", "total_bytes"])
    return df


def parse_count_abun_table(count_abun_table):
    sample_id = os.path.basename(count_abun_table).split(".")[0]

    df = pd.read_csv(count_abun_table, sep="\t").set_index("CAZy_gene_name")
    df_count = df.loc[:, ["CAZy_gene_copy_number"]]\
                 .rename(columns={"CAZy_gene_copy_number": sample_id})

    df_abun = df.loc[:, ["CAZy_gene_relab_abun"]]\
                .rename(columns={"CAZy_gene_relab_abun": sample_id})
    return df_count, df_abun


def merger(args):
    count_abun_list = [line.strip() for line in xopen(args.tables).readlines()]
    df_count_list = []
    df_abun_list = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
       for df_count, df_abun in executor.map(parse_count_abun_table, count_abun_list):
           df_count_list.append(df_count)
           df_abun_list.append(df_abun)

    count_df = pd.concat(df_count_list, axis=1, sort=True).fillna(0).reset_index()
    abun_df = pd.concat(df_abun_list, axis=1, sort=True).fillna(0).reset_index()

    count_df.to_csv(args.output_count_profile, sep="\t", index=False)
    abun_df.to_csv(args.output_abun_profile, sep="\t", index=False)
    return count_df, abun_df 


def profiler(args):
    cazy_dict = {}
    record_dict = {}
    first_seq_id = ""
    common_list = []

    with xopen(args.tables) as handle:
        blast_record = handle.readline().strip().split("\t")
        seq_id, pair = blast_record[0].split("/")

        record_dict[pair] = set() 
        record_dict[pair].add((blast_record[1], blast_record[2], blast_record[10]))

        first_seq_id = seq_id

        for line in handle:
            blast_record = line.strip().split("\t")
            seq_id, pair = blast_record[0].split("/")
            if seq_id == first_seq_id:
                if pair in record_dict:
                    record_dict[str(pair)].add((blast_record[1], float(blast_record[2]), float(blast_record[10])))
                else:
                    record_dict[str(pair)] = set() 
                    record_dict[str(pair)].add((blast_record[1], float(blast_record[2]), float(blast_record[10])))
            else:
                if len(record_dict) == 2:
                    common_ref = record_dict["1"] & record_dict["2"]
                    common_ref_len = len(common_ref)
                    if common_ref_len != 0:
                        sseqid = ""
                        pident = 0.0
                        evalue = 100
                        for sseq in common_ref:
                            if (sseq[1] >= pident) & (sseq[2] <= evalue):
                                sseqid = sseq[0]
                                pident = sseq[1]
                                evalue = sseq[2]
                        if sseqid in cazy_dict:    
                            cazy_dict[sseqid] += 1
                        else:
                            cazy_dict[sseqid] = 1
                        common_list.append([first_seq_id, common_ref_len, sseqid, pident, evalue])

                record_dict = {}
                first_seq_id = seq_id
                record_dict[str(pair)] = set() 
                record_dict[str(pair)].add((blast_record[1], blast_record[2], blast_record[10]))
 
        if len(record_dict) == 2:
            common_ref = record_dict["1"] & record_dict["2"]
            common_ref_len = len(common_ref)
            if common_ref_len != 0:
                sseqid = ""
                pident = 0
                evalue = 1000000
                for sseq in common_ref:
                    if (sseq[1] >= pident) and (sseq[2] <= evalue):
                        sseqid = sseq[0]
                        pident = sseq[1]
                        evalue = sseq[2]
                if sseqid in cazy_dict:    
                    cazy_dict[sseqid] += 1
                else:
                    cazy_dict[sseqid] = 1
                common_list.append([first_seq_id, common_ref_len, sseqid, pident, evalue])


    mapping_count_table = pd.DataFrame(common_list,
                                       columns=["reads_name", "r1_r2_common_ref_number",
                                                "ref_name", "identity", "evalue"])
    mapping_count_table.to_csv(args.mapping_count_table, sep="\t", index=False)

    cazy_fai_df = parse_fai(args.cazydb_fai)
    #print(cazy_fai_df.head())
    
    #print(cazy_dict.items())
    cazy_df = pd.DataFrame(cazy_dict.items(), columns=["CAZy_gene_name", "CAZy_gene_copy_number"])\
        .merge(cazy_fai_df).loc[:, ["CAZy_gene_name", "CAZy_gene_len", "CAZy_gene_copy_number"]]
    #print(cazy_df)
    
    cazy_df["CAZy_gene_copy_number_norm"] = cazy_df["CAZy_gene_copy_number"] / cazy_df["CAZy_gene_len"]
    cazy_df["CAZy_gene_relab_abun"] = cazy_df["CAZy_gene_copy_number_norm"] / sum(cazy_df["CAZy_gene_copy_number_norm"])
    cazy_df.to_csv(args.cazy_abundance_table, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
       prog="CAZyme_profiler.py",
       usage="CAZyme_profiler.py [subcommand] [options]",
       description="CAZyme gene profiler based on mapping metagenomcis reads to CAZydb using diaomond blastx"
    )

    parser.add_argument(
        "-v",
        "--version",
        default=False,
        help="print version and exit"
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "--tables",
        default=sys.stdin,
        help="files to read, blast-m6-table or profile-table, if empty, stdin will be used"
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")

    parser_profiler = subparsers.add_parser(
        "profiler",
        parents=[parent_parser],
        prog="CAZyme_profiler.py profiler",
        description="blastout parser, calculate gene copy number and relative abundance",
        help="parse blastout to calculate CAZy gene profile"
    )
    parser_profiler.add_argument("--cazydb-fai", dest="cazydb_fai", help="CAZydb length info")
    parser_profiler.add_argument("--mapping-count-table", dest="mapping_count_table", help="reads pair mapping count table")
    parser_profiler.add_argument("--cazy-abundance-table", dest="cazy_abundance_table", help="CAZy abundance table")

    parser_merger = subparsers.add_parser(
        "merger",
        parents=[parent_parser],
        prog="CAZyme_profiler.py merger",
        help="merge CAZyme count and abundance tables"
    )
    parser_merger.add_argument(
        '--output-count-profile',
        dest="output_count_profile",
        help="output_count_profile"
    )
    parser_merger.add_argument(
        '--output-abun-profile',
        dest="output_abun_profile",
        help="output_abun_profile"
    )
    parser_merger.add_argument(
        "--threads", type=int, default=8, help="threads, default: 8"
    )

    parser_profiler._optionals.title = "arguments"
    parser_profiler.set_defaults(func=profiler)
    parser_merger._optionals.title = "arguments"
    parser_merger.set_defaults(func=merger)

    args = parser.parse_args()
    try:
        if args.version:
            print("CAZyme_profiler.py v0.1.0")
            sys.exit(0)
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == '__main__':
    main()