#!/usr/bin/env python3

import os
from Bio import SeqIO

import pandas as pd


def get_binning_info(bins_dir, cluster_file, assembler):
    if assembler.lower() in ["spades", "metaspades", "megahit"]:
        with os.scandir(bins_dir) as itr, open(cluster_file, "w") as oh:
            for entry in itr:
                bin_id, suffix = os.path.splitext(entry.name)
                if suffix == ".fa":
                    cluster_num = bin_id.split(".")[-1]
                    bin_fa = os.path.join(bins_dir, entry.name)
                    for seq in SeqIO.parse(bin_fa, "fasta"):
                        # graphbin
                        # oh.write("%s,%s" %
                        #         ("_".join(seq.id.split("_")[:2]), cluster_num))
                        # graphbin 2
                        oh.write(f"{seq.id},{cluster_num}\n")


def generate_bins(cluster_file, scaftigs, prefix):

    def get_accession(identifier):
        return "_".join(identifier.split("_")[:2])

    # graphbin
    # scaftigs_index = SeqIO.index(scaftigs, "fasta", key_function=get_accession)

    # graphbin2
    scaftigs_index = SeqIO.index(scaftigs, "fasta")

    df = pd.read_csv(cluster_file, names=["scaftigs_id", "bin_id"])\
           .astype({"scaftigs_id": str,
                    "bin_id": str})\
           .set_index("bin_id")

    for i in df.index.unique():
        scaftigs_id_list = df.loc[[i], "scaftigs_id"]\
                             .dropna().tolist()
        bin_fa = prefix + "." + i + ".fa"
        with open(bin_fa, 'w') as oh:
            for scaftigs_id in scaftigs_id_list:
                SeqIO.write(scaftigs_index[scaftigs_id], oh, "fasta")


def extract_bins_report(bins_report_table):
    bins_report = pd.read_csv(bins_report_table, sep='\t', header=[0, 1])\
                    .rename(columns={
                            "Unnamed: 0_level_1": "assembly_group",
                            "Unnamed: 1_level_1": "bin_id",
                            "Unnamed: 2_level_1": "bin_file",
                            "Unnamed: 3_level_1": "assembler",
                            "Unnamed: 4_level_1": "binner"}, level=1)

    bins_report = bins_report[[
        ("assembly_group", "assembly_group"),
        ("bin_id", "bin_id"),
        ("bin_file", "bin_file"),
        ("assembler", "assembler"),
        ("binner", "binner"),
        ("length", "sum"),
        ("length", "N50")]]
    
    bins_report.columns = ["assembly_group", "bin_id", "bin_file", "assembler", "binner", "length", "N50"]
    return bins_report


'''
            table_bins = pd.read_csv(input.table_bins, sep="\t", header=[0, 1])
            table_bins = table_bins[
                [
                    ("bin_id", "Unnamed: 1_level_1"),
                    ("chr", "count"),
                    ("length", "sum"),
                    ("length", "min"),
                    ("length", "max"),
                    ("length", "std"),
                    ("length", "N50")
                ]
            ]
            table_bins.columns = [
                "user_genome",
                "contig_number",
                "contig_length_sum",
                "contig_length_min",
                "contig_length_max",
                "contig_length_std",
                "N50"
            ]
'''