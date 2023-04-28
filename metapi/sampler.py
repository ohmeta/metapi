#!/usr/bin/env python
import glob
import os
import sys

import pandas as pd


def pass_gzip_format(fq_list):
    cancel = False
    for fq in fq_list:
        if pd.isna(fq):
            pass
        elif not fq.endswith(".gz"):
            cancel = False
            print(f"{fq} is not gzip format")
    return cancel


def parse_samples(samples_tsv):
    samples_df = pd.read_csv(samples_tsv, sep="\t")
    samples_df = samples_df.applymap(str)

    samples_df["multibinning_group"] = samples_df["assembly_group"] + "__" + samples_df["binning_group"]
    samples_df = samples_df.set_index(["sample_id", "assembly_group", "binning_group", "multibinning_group"])

    cancel = False

    sra_headers = [
        "sra_pe",
        "sra_se",
        "sra_long"
    ]

    fastq_headers = [
        "short_forward_reads",
        "short_reverse_reads",
        "short_intervead_reads",
        "short_single_reads",
        "long_reads"
    ]

    # check header
    is_sra = False
    is_fastq = False
    for sra in sra_headers:
        if sra in samples_df.columns:
            is_sra = True
    for fq in fastq_headers:
        if fq in samples_df.columns:
            is_fastq = True
    if is_sra and is_fastq:
        cancel = True
        print("Can't process sra and fastq at the same time, please provide fastq only or sra only")
    if cancel:
        sys.exit(-1)

    # check sample_id
    for sample_id in samples_df.index.get_level_values("sample_id").unique():
        if "." in sample_id:
            if "." in sample_id:
                print(f"{sample_id} contain '.', please remove '.'")
                cancel = True
    if cancel:
        print(f"{samples_tsv} sample_id contain '.' please remove '.', now quiting :)")

    # check fastq headers
    if is_fastq:
        if ("short_forward_reads" in samples_df.columns) and ("short_reverse_reads" in samples_df.columns):
            if "short_intervedad_reads" in samples_df.columns:
                cancel = True
                print("can't specific short_forward_reads, short_reverse_reads and short_intervead_reads at the same time")
                print("please only provide short_forward_reads and short_reverse_reads, or short_intervead_reads")
            else:
                pass
        elif ("short_forward_reads" in samples_df.columns) or ("short_rverse_reads" in samples_df.columns):
            cancel = True
            print("Please specific short_forward_reads and short_reverse_reads at the same time")
        else:
            pass

    # check reads_format
    for sample_id in samples_df.index.get_level_values("sample_id").unique():
        # check short paired-end reads
        if ("short_forward_reads" in samples_df.columns) and ("short_reverse_reads" in samples_df.columns):
            forward_reads = samples_df.loc[sample_id, :, :, :]["short_forward_reads"].tolist()
            reverse_reads = samples_df.loc[sample_id, :, :, :]["short_reverse_reads"].tolist()

            for forward_read, reverse_read in zip(forward_reads, reverse_reads):
                if pd.isna(forward_read) and pd.isna(reverse_read):
                    pass
                elif (not pd.isna(forward_read)) and (not pd.isna(reverse_read)):
                    if forward_read.endswith(".gz") and reverse_read.endswith(".gz"):
                        pass
                    else:
                        cancel = True
                        print(f"{forward_read} or {reverse_read} is not gzip format")
                else:
                    cancel = True
                    print(f"It seems short paired-end reads only specific forward or reverse reads, please check again!")

        # check short single-end reads
        if "short_single_reads" in samples_df.columns:
            single_reads = samples_df.loc[sample_id, :, :, :]["short_single_reads"].tolist()
            cancel = pass_gzip_format(single_reads)

        # check long reads
        if "long_reads" in samples_df.columns:
            long_reads = samples_df.loc[sample_id, :, :, :]["long_reads"].tolist()
            cancel = pass_gzip_format(long_reads)

        # check short_interleaved
        if "short_interleaved" in samples_df.columns:
            interleaved_reads = samples_df.loc[sample_id, :, :, :]["long_reads"].tolist()
            cancel = pass_gzip_format(interleaved_reads)

    if cancel:
        sys.exit(-1)
    else:
        return samples_df


def parse_mags(mags_dir):
    bin_list = []
    for bin_ in glob.glob(mags_dir + "/*/*.fa"):
        bin_dict = dict()
        bin_dict["path"] = bin_.strip()
        bin_dict["id"] = os.path.basename(bin_).rstrip(".fa")
        bin_list.append(bin_dict)
    mags_df = pd.DataFrame(bin_list).set_index("id", drop=False)
    return mags_df


def get_reads(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, :, :, :][col].dropna().tolist()


def get_samples_id_by_assembly_group(sample_df, assembly_group):
    return sample_df.loc[:, assembly_group, :, :].index.get_level_values("sample_id").unique()


def get_samples_id_by_binning_group(sample_df, binning_group):
    return sample_df.loc[:, :, binning_group, :].index.get_level_values("sample_id").unique()


def get_samples_id_by_assembly_and_binning_group(sample_df, assembly_group, binning_group):
    return sample_df.loc[:, assembly_group, binning_group, :].index.get_level_values("sample_id").unique()


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