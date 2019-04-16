#!/usr/bin/env python
import glob
import os
import pandas as pd


def samples_validator(sample_df, input_type):
    error_count = 0
    for i in sample_df.index:
        if input_type == "fastq":
            fq1 = sample_df.loc[[i], "fq1"].dropna().tolist()
            fq2 = sample_df.loc[[i], "fq2"].dropna().tolist()
            for r1, r2 in zip(fq1, fq2):
                if (not os.path.exists(r1)) or (not os.path.exists(r2)):
                    print("error:\t%s\t%s\t%s" % (i, r1, r2))
                    error_count += 1
        elif input_type == "sra":
            for sra in sample_df.loc[[i], "sra"].dropna().tolist():
                if not os.path.exists(sra):
                    print("error:\t%s\t%s" % (i, sra))
                    error_count += 1
        else:
            print("wrong input type! just support fastq or sra")
    return error_count


def parse_samples(samples_tsv, input_type, check=True):
    samples_df = pd.read_csv(samples_tsv, sep='\s+').set_index("id", drop=False)
    if check:
        error_count = samples_validator(samples_df, input_type)
        if error_count == 0:
            return samples_df
        else:
            print("find %d error" % error_count)
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


def parse_cobin_samples_id(query_list):
    with open(query_list, 'r') as ih:
        samples_id = [line.strip() for line in ih]
    return samples_id


def renamed_id(samples_df, wildcards):
    return samples_df.loc[[wildcards.sample], "id_2"].dropna().tolist()[0]
