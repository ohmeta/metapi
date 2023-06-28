#!/usr/bin/env python3

import time
import gzip
import os
import sys
import subprocess

import pandas as pd
from Bio import bgzf


def gtdbtk_prepare_from_mags(rep_table, batch_num, mags_dir):
    os.makedirs(mags_dir, exist_ok=True)

    table_df = pd.read_csv(rep_table, sep="\t").query('MIMAG_quality_level!="low_quality"').reset_index(drop=True)

    batchid = -1
    if len(table_df) > 0:
        for batch in range(0, len(table_df), batch_num):
            batchid += 1
            columns_list = list(table_df.columns)

            bin_file_index = columns_list.index("bin_file")
            genome_index = columns_list.index("genome")
            best_translation_table_index = columns_list.index("best_translation_table")

            table_split = table_df.iloc[batch:batch+batch_num,
                                        [bin_file_index, genome_index, best_translation_table_index]]

            table_split\
                .to_csv(os.path.join(mags_dir, f"mags_input.{batchid}.tsv"),
                        sep="\t", index=False, header=None)
    else:
        subprocess.run(f'''touch {os.path.join(mags_dir, "mags_input.0.tsv")}''', shell=True)


def gtdbtk_prepare_from_genes(rep_table, batch_num, mags_dir):
    os.makedirs(mags_dir, exist_ok=True)

    table_df = pd.read_csv(rep_table, sep="\t")

    batchid = -1
    if len(table_df) > 0:
        for batch in range(0, len(table_df), batch_num):
            batchid += 1
            columns_list = list(table_df.columns)

            pep_file_index = columns_list.index("pep_file")
            best_translation_table_index = columns_list.index("best_translation_table")

            table_split = table_df.iloc[batch:batch+batch_num,
                                        [pep_file_index, best_translation_table_index]]

            table_split["pep_location"] = table_split.apply(lambda x: os.path.splitext(x["pep_file"])[0], axis=1)
            table_split["pep_basename"] = table_split.apply(lambda x: os.path.basename(x["pep_location"]), axis=1)

            table_split\
                .loc[:, ["pep_location", "pep_basename", "best_translation_table"]]\
                .to_csv(os.path.join(mags_dir, f"mags_input.{batchid}.tsv"),
                        sep="\t", index=False, header=None)
            pepcount = 0 
            for pep_file in table_split["pep_file"]:
                pep_file_ = os.path.splitext(pep_file)[0]
                if not os.path.exists(pep_file_):
                    subprocess.run(f"pigz -dkf {pep_file}", shell=True)
                if os.path.exists(pep_file_):
                    pepcount += 1
            if pepcount == len(table_split):
                print(f"Uncompress done for batchid: {batchid}")
            else:
                print(f"Uncompress failed for batchid: {batchid}")
                print(f"Please check when prepare input for gtdbtk")
                sys.exit(-1)
    else:
        subprocess.run(f'''touch {os.path.join(mags_dir, "mags_input.0.tsv")}''', shell=True)


def demultiplex(kraken2_output, r1, r2, change_seq_id, prefix, log=None):
    start_time = time.time()
    taxid_counter = {}
    demultiplexer = {}

    with open(kraken2_output, "r") as kh:
        for line in kh:
            cols = line.split("\t")
            read_id = cols[1]
            tax_name = cols[2].split("(")[0].strip()
            tax_id = int(cols[2].split("(")[-1].split(")")[0].split()[-1])

            demultiplexer[read_id] = tax_id
            if tax_id in taxid_counter:
                taxid_counter[tax_id][1] += 1
            else:
                taxid_counter[tax_id] = [tax_name, 1]
    if log is not None:
        log_h.write(
            "step_1: parse kraken2 output has spent %d s\n" % (time.time() - start_time)
        )

    start_time = time.time()
    gzip_h = {}
    for i in taxid_counter:
        gzip_h[i] = {}
        gzip_h[i]["r1"] = bgzf.BgzfWriter(prefix + ".%d.1.fq.gz" % i, "wb")
        gzip_h[i]["r2"] = bgzf.BgzfWriter(prefix + ".%d.2.fq.gz" % i, "wb")

    if r1.endswith(".gz"):
        r1_h = gzip.open(r1, "rt")
        r2_h = gzip.open(r2, "rt")
    else:
        r1_h = open(r1, "rt")
        r2_h = open(r2, "rt")

    if change_seq_id:
        sample_tag = os.path.basename(prefix)

    if log is not None:
        log_h.write("step_2: begin demultiplex taxid-reads\n")
    for read_1, read_2 in zip(r1_h, r2_h):
        read_id = read_1[1:].split("/")[0]
        if change_seq_id:
            gzip_h[demultiplexer[read_id]]["r1"].write(
                ">%s|%s%s%s%s"
                % (sample_tag, read_1[1:], next(r1_h), next(r1_h), next(r1_h))
            )
            gzip_h[demultiplexer[read_id]]["r2"].write(
                ">%s|%s%s%s%s"
                % (sample_tag, read_2[1:], next(r2_h), next(r2_h), next(r2_h))
            )
        else:
            gzip_h[demultiplexer[read_id]]["r1"].write(
                "%s%s%s%s" % (read_1, next(r1_h), next(r1_h), next(r1_h))
            )
            gzip_h[demultiplexer[read_id]]["r2"].write(
                "%s|%s%s%s%s" % (read_2, next(r2_h), next(r2_h), next(r2_h))
            )
    if log is not None:
        log_h.write(
            "step_2: demultiplex taxid-reads has spent %d s\n"
            % (time.time() - start_time)
        )
