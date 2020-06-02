#!/usr/bin/env python
import glob
import gzip
import os
import sys
import time
import pandas as pd
from Bio import bgzf


def parse_samples(config):
    samples_df = pd.read_csv(config["params"]["samples"], sep="\s+").set_index(
        "id", drop=False
    )

    cancel = False
    if "fq1" in samples_df.columns:
        for sample_id in samples_df.index.unique():
            fq1_list = samples_df.loc[[sample_id], "fq1"].dropna().tolist()
            for fq_file in fq1_list:
                if not fq_file.endswith(".gz"):
                    print("%s need gzip format" % fq_file)
                    cancel = True
                if not os.path.exists(fq_file):
                    print("%s not exists" % fq_file)
                    cancel = True

    if cancel:
        sys.exit(0)
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
