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
