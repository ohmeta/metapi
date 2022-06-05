#!/usr/bin/env python
import glob
import gzip
import os
import sys
import time

import pandas as pd
from Bio import bgzf


def parse_samples(
    samples_tsv,
    interleaved="false",
    reads_layout="pe",
    begin_point="trimming",
    check_samples=False,
):
    samples_df = pd.read_csv(
        samples_tsv, sep="\t", dtype={
            "sample_id": str,
            "fq1": str, "fq2": str,
            "assembly_group": str,
            "binning_group": str})
    samples_df["multibinning_group"] = samples_df["assembly_group"] + "__" + samples_df["binning_group"]
    samples_df = samples_df.set_index(["sample_id", "assembly_group", "binning_group", "multibinning_group"])

    cancel = False
    if "fq1" in samples_df.columns:
        for sample_id in samples_df.index.get_level_values("sample_id").unique():
            sample_id = str(sample_id)
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.', now quiting :)")
                cancel = True

            fq1_list = samples_df.loc[sample_id, :, :, :]["fq1"].dropna().tolist()
            fq2_list = samples_df.loc[sample_id, :, :, :]["fq2"].dropna().tolist()
            for fq_file in fq1_list:
                if not fq_file.endswith(".gz"):
                    print(f"{fq_file} need gzip format")
                    cancel = True
                if check_samples:
                    if not os.path.exists(fq_file):
                        print(f"{fq_file} not exists")
                        cancel = True
                    if (reads_layout == "pe") and (not interleaved):
                        if len(fq2_list) == 0:
                            print(f"{sample_id} fq2 not exists")
                            cancel = True
    elif "sra" in samples_df.columns:
        for sample_id in samples_df.index.get_level_values("sample_id").unique():
            sample_id = str(sample_id)
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.', now quiting :)")
                cancel = True

            if check_samples:
                sra_list = samples_df.loc[sample_id, :, :, :]["sra"].dropna().tolist()
                for sra_file in sra_list:
                    if not os.path.exists(sra_file):
                        print(f"{sra_file} not exists")
                        cancel = True
    else:
        print("wrong header: {header}".format(header=samples_df.columns))
        cancel = True

    if begin_point == "binning":
        if len(samples_df) != len(samples_df.index.unique()):
            print(
                "when begin with binning, samples id need to be unique, because we can't merge assembly"
            )
            cancel = True

        if check_samples:
            if "scaftigs" in samples_df.columns:
                for sample_id in samples_df.index.unique():
                    sample_id = str(sample_id)
                    scaftigs = samples_df.loc[sample_id, "scaftigs"]
                    if not os.path.exists(scaftigs):
                        print(f"{scaftigs} not exists")
                        cancel = True

    if cancel:
        sys.exit(-1)
    else:
        return samples_df


def parse_bins(bins_dir):
    bin_list = []
    for bin_ in glob.glob(bins_dir + "/*/*.fa"):
        bin_dict = dict()
        bin_dict["path"] = bin_.strip()
        bin_dict["id"] = os.path.basename(bin_).rstrip(".fa")
        bin_list.append(bin_dict)
    bins = pd.DataFrame(bin_list).set_index("id", drop=False)
    return bins


def get_reads(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, :, :, :][col].dropna().tolist()


def get_samples_id_by_assembly_group(sample_df, assembly_group):
    return sample_df.loc[:, assembly_group, :, :].index.get_level_values("sample_id").unique()


def get_samples_id_by_binning_group(sample_df, binning_group):
    return sample_df.loc[:, :, binning_group, :].index.get_level_values("sample_id").unique()


def get_assembly_group_by_binning_group(sample_df, binning_group):
    return sample_df.loc[:, :, binning_group, :].index.get_level_values("assembly_group").unique()


def get_binning_group_by_assembly_group(sample_df, assembly_group):
    return sample_df.loc[:, assembly_group, :, :].index.get_level_values("binning_group").unique()


def get_multibinning_group_by_assembly_group(sample_df, assembly_group):
    return sample_df.loc[:, assembly_group, :, :].index.get_level_values("multibinning_group").unique()


def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, :, :][col].dropna()[0]


def get_sample_id_(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample_, :, :][col].dropna()[0]


def get_bin_id(bin_df, wildcards, col):
    return bin_df.loc[wildcards.bin, [col]].dropna()[0]


def sample_cluster_info(sample_df):
    '''
    vamb cluster
    '''
    binning_groups = sample_df.index.get_level_values("binning_group").unique()

    for binning_group in binning_groups:
        assembly_groups = sample_df.loc[:, :, binning_group].index.get_level_values("assembly_group").unique()
        assembly_groups = sorted(assembly_groups)
        for i in range(0, len(assembly_groups)):
            print(f'''{assembly_groups[i]}\tS{i+1}''')