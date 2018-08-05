#!/usr/bin/env python

import glob
import os

import pandas


def parse_samples(samples_tsv):
    samples = pandas.read_table(samples_tsv).set_index("id", drop=False)
    return samples


def parse_bins(bins_dir):
    bin_list = []
    for bin in glob.glob(bins_dir + "/*/*bin*fa"):
        bin_dict = {}
        bin_dict["path"] = bin.strip()
        bin_dict["id"] = os.path.basename(bin).rstrip(".fa")
        bin_list.append(bin_dict)
    bins = pandas.DataFrame(bin_list).set_index("id", drop=False)
    return bins


def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]


def get_bin_id(bin_df, wildcards, col):
    return bin_df.loc[wildcards.bin, [col]].dropna()[0]
