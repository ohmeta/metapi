#!/usr/bin/env python
import glob
import gzip
import os
import sys
import time
import pandas as pd
from Bio import bgzf


def check_duplicated(samples_df, begin, format):
    if begin == "assembly":
        cancel_task = False
        if format == "fastq":
            for i in samples_df.index.unique():
                count = len(samples_df.loc[[i], "fq1"].dropna().tolist())
                if count > 1:
                    cancel_task = True
                    print("exists duplicated sample: %s" % i)
        if cancel_task:
            sys.exit(1)


def check_exists(samples_df, input_type, is_pe):
    error_count = 0
    for i in samples_df.index:
        if input_type == "fastq":
            if is_pe:
                fq1 = samples_df.loc[[i], "fq1"].dropna().tolist()
                fq2 = samples_df.loc[[i], "fq2"].dropna().tolist()
                for r1, r2 in zip(fq1, fq2):
                    if (not os.path.exists(r1)) or (not os.path.exists(r2)):
                        print("Error:\t%s\t%s\t%s" % (i, r1, r2))
                        error_count += 1
            else:
                fq = samples_df.loc[[i], "fq1"].dropna().tolist()
                for r in fq:
                    if not os.path.exists(r):
                        print("Error:\t%s\t%s" % (i, r))
                        error_count += 1
        elif input_type == "sra":
            for sra in samples_df.loc[[i], "sra"].dropna().tolist():
                if not os.path.exists(sra):
                    print("Error:\t%s\t%s" % (i, sra))
                    error_count += 1
        else:
            print("Wrong input type! just support fastq or sra")
    return error_count


def parse_samples(config):
    samples_df = pd.read_csv(config["params"]["samples"], sep="\s+").set_index(
        "id", drop=False
    )
    is_pe = True if config["params"]["reads_layout"] == "pe" else False
    error_count = check_exists(samples_df, config["params"]["reads_format"], is_pe)
    if error_count == 0:
        check_duplicated(samples_df, config["params"]["begin"], format)
        return samples_df
    else:
        print("find %d error" % error_count)
        sys.exit(1)
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
    print(wildcards.sample)
    return sample_df.loc[[wildcards.sample], col].dropna().tolist()


def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]


def get_sample_id_(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample_, [col]].dropna()[0]


def get_bin_id(bin_df, wildcards, col):
    return bin_df.loc[wildcards.bin, [col]].dropna()[0]


def parse_cobin_samples_id(samples_df, config):
    samples_id = []
    if config["params"]["cobinning"]["do"]:
        with open(config["params"]["cobinning"]["renamed_id"], "r") as ih:
            samples_id = [line.strip() for line in ih]
        if config["params"]["cobinning"]["rename"]:
            rename_dict = {
                k: "S" + str(v + 1)
                for k, v in zip(
                    samples_df.index.unique(), range(len(samples_df.index.unique()))
                )
            }
            samples_df = samples_df.assign(
                id_2=samples_df.id.apply(lambda x: rename_dict[x])
            )
            samples_df.loc[:, ["id", "id_2"]].drop_duplicates(["id", "id_2"]).to_csv(
                config["params"]["cobinning"]["rename_id"], sep="\t", index=False
            )
    return samples_id


def renamed_id(samples_df, wildcards):
    return samples_df.loc[[wildcards.sample], "id_2"].dropna().tolist()[0]
