#!/usr/bin/env python
import glob
import gzip
import os
import sys
import time

import pandas as pd
from Bio import bgzf


def parse_samples(config):
    samples_df = pd.read_csv(config["params"]["samples"], sep="\t").set_index(
        "id", drop=False
    )

    cancel = False
    if "fq1" in samples_df.columns:
        for sample_id in samples_df.index.unique():
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.', now quiting :)")
                cancel = True
            fq1_list = samples_df.loc[[sample_id], "fq1"].dropna().tolist()
            fq2_list = samples_df.loc[[sample_id], "fq2"].dropna().tolist()
            for fq_file in fq1_list:
                if not fq_file.endswith(".gz"):
                    print(f"{fq_file} need gzip format")
                    cancel = True
                if not os.path.exists(fq_file):
                    print(f"{fq_file} not exists")
                    cancel = True
                if (config["params"]["reads_layout"] == "pe") and (
                    not config["params"]["interleaved"]
                ):
                    if len(fq2_list) == 0:
                        print(f"{sample_id} fq2 not exists")
                        cancel = True
    elif "sra" in samples_df.columns:
        for sample_id in samples_df.index.unique():
            sra_list = samples_df.loc[[sample_id], "sra"].dropna().tolist()
            for sra_file in sra_list:
                if not os.path.exists(sra_file):
                    print(f"{sra_file} not exists")
                    cancel = True
    else:
        print("wrong header: {header}".format(header=samples_df.columns))
        cancel = True

    if cancel:
        sys.exit(-1)
    else:
        return samples_df


def parse_bins(bins_dir):
    bin_list = []
    for bin_ in glob.glob(bins_dir + "/*/*bin*fa"):
        bin_dict = dict()
        bin_dict["path"] = bin_.strip()
        bin_dict["id"] = os.path.basename(bin_).rstrip(".fa")
        bin_list.append(bin_dict)
    bins = pd.DataFrame(bin_list).set_index("id", drop=False)
    return bins


def get_reads(sample_df, wildcards, col):
    return sample_df.loc[[wildcards.sample], col].dropna().tolist()


def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]


def get_sample_id_(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample_, [col]].dropna()[0]


def get_bin_id(bin_df, wildcards, col):
    return bin_df.loc[wildcards.bin, [col]].dropna()[0]
