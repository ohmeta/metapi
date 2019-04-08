#!/usr/bin/env python

import glob
import os
import sys

import pandas


def samples_validator(sample_df):
    error_count = 0
    for i in sample_df.index:
        fq1, fq2 = sample_df.loc[i, ["fq1", "fq2"]]
        if (not os.path.exists(fq1)) or (not os.path.exists(fq2)):
            print("error:\t%s\t%s\t%s" % (i, fq1, fq2))
            error_count += 1
    return error_count


def parse_samples(samples_tsv, check=True):
    samples_df = pandas.read_csv(samples_tsv, sep='\s+').set_index("id", drop=False)
    if check:
        error_count = samples_validator(samples_df)
        if error_count == 0:
            return samples_df
        else:
            print("find %d error" % error_count)
    else:
        return samples_df


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


def get_sample_id_(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample_, [col]].dropna()[0]


def get_bin_id(bin_df, wildcards, col):
    return bin_df.loc[wildcards.bin, [col]].dropna()[0]


def parse_cobin_samples_id(query_list):
    samples_id = []
    with open(query_list, 'r') as ih:
        samples_id = [line.strip() for line in ih]
    return samples_id
