#!/usr/bin/env python3

import os
from Bio import SeqIO

import pandas as pd


def get_binning_info(bins_dir, cluster_file, bin_suffix, assembler):
    if assembler.lower() in ["spades", "metaspades", "megahit"]:
        with os.scandir(bins_dir) as itr, open(cluster_file, "w") as oh:
            for entry in itr:
                bin_id, suffix = os.path.splitext(entry.name)
                if suffix == "." + bin_suffix:
                    cluster_num = bin_id.split(".")[-1]
                    bin_fa = os.path.join(bins_dir, entry.name)
                    for seq in SeqIO.parse(bin_fa, "fasta"):
                        # graphbin
                        # oh.write("%s,%s" %
                        #         ("_".join(seq.id.split("_")[:2]), cluster_num))
                        # graphbin 2
                        oh.write(f"{seq.id},{cluster_num}\n")


def generate_bins(cluster_file, scaftigs, prefix, suffix):

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
        bin_fa = prefix + "." + i + "." + suffix
        with open(bin_fa, 'w') as oh:
            for scaftigs_id in scaftigs_id_list:
                SeqIO.write(scaftigs_index[scaftigs_id], oh, "fasta")
