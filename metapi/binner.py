#!/usr/bin/env python3

import os
from Bio import SeqIO

import pandas as pd


def get_binning_info(bins_dir, cluster_file, bin_suffix):
    cluster_dict = {}

    with os.scandir(bins_dir) as itr, open(cluster_file, "w") as oh:
        for entry in itr:
            bin_id, bin_suffix = os.path.splitext(entry.name)
            if bin_suffix == "." + bin_suffix:
                cluster_num = bin_id.split(".")[-1]
                for seq in SeqIO.parse(os.path.join(bins_dir, entry.name), "fasta"):
                    oh.write("%s,%s\n" % (seq.id, bin_id))


def generate_bins(cluster_file, scaftigs, prefix, suffix):
    scaftigs_index = SeqIO.index(scaftigs, "fasta")
    df = pd.read_csv(cluster_file, names=["scaftigs_id", "cluster_id"])\
           .set_index("cluster_id")\
           .as_type({"scaftigs_id": str,
                     "cluster_id": str})

    for i in df.index.unique():
        scaftigs_id_list = df.loc[["cluster_id"], "scaftigs_id"].dropna().tolist()
        bin_fa = prefix + "." + i + "." + suffix
        with open(bin_fa, 'w') as oh:
            SeqIO.write(scaftigs_index[i], oh, "fasta")
