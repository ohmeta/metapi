#!/usr/bin/env python
import glob
import os
import sys

import pandas as pd


SRA_HEADERS = {
    "PE": "sra_pe",
    "SE": "sra_se",
    "LONG": "sra_long"
}

FQ_HEADERS = {
    "PE_FORWARD": "short_forward_reads",
    "PE_REVERSE": "short_reverse_reads",
    "INTERLEAVED": "short_interleaved_reads",
    "SE": "short_single_reads",
    "LONG": "long_reads"
}

HEADERS = {
    "SRA": SRA_HEADERS,
    "FQ": FQ_HEADERS
}


def pass_gzip_format(fq_list):
    cancel = False
    for fq in fq_list:
        if pd.isna(fq):
            pass
        elif not fq.endswith(".gz"):
            cancel = True
            print(f"{fq} is not gzip format")
    return cancel


def parse_samples(samples_tsv):
    samples_df = pd.read_csv(samples_tsv, sep="\t")
    #samples_df = samples_df.applymap(str)

    samples_df["multibinning_group"] = samples_df["assembly_group"] + "__" + samples_df["binning_group"]
    samples_df = samples_df.set_index(["sample_id", "assembly_group", "binning_group", "multibinning_group"])

    cancel = False

    # check header
    is_sra = False
    is_fastq = False
    for sra_k, sra_v in SRA_HEADERS.items():
        if sra_v in samples_df.columns:
            is_sra = True
    for fq_k, fq_v in FQ_HEADERS.items():
        if fq_v in samples_df.columns:
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
        if (FQ_HEADERS["PE_FORWARD"] in samples_df.columns) and (FQ_HEADERS["PE_REVERSE"] in samples_df.columns):
            if FQ_HEADERS["INTERLEAVED"] in samples_df.columns:
                cancel = True
                print(f'''can't specific {FQ_HEADERS["PE_FORWARD"]}, {FQ_HEADERS["PE_REVERSE"]} and {FQ_HEADERS["INTERLEAVED"]} at the same time''')
                print(f'''please only specific {FQ_HEADERS["PE_FORWARD"]} and {FQ_HEADERS["PE_REVERSE"]}, or {FQ_HEADERS["INTERLEAVED"]}''')
            else:
                pass
        elif (FQ_HEADERS["PE_FORWARD"] in samples_df.columns) or (FQ_HEADERS["PE_REVERSE"] in samples_df.columns):
            cancel = True
            print(f'''please only specific {FQ_HEADERS["PE_FORWARD"]} and {FQ_HEADERS["PE_REVERSE"]} at the same time''')
        else:
            pass

    # check reads_format
    for sample_id in samples_df.index.get_level_values("sample_id").unique():
        # check short paired-end reads
        if (FQ_HEADERS["PE_FORWARD"] in samples_df.columns) and (FQ_HEADERS["PE_REVERSE"] in samples_df.columns):
            forward_reads = samples_df.loc[sample_id, :, :, :][FQ_HEADERS["PE_FORWARD"]].tolist()
            reverse_reads = samples_df.loc[sample_id, :, :, :][FQ_HEADERS["PE_REVERSE"]].tolist()

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
        if FQ_HEADERS["SE"] in samples_df.columns:
            single_reads = samples_df.loc[sample_id, :, :, :][FQ_HEADERS["SE"]].tolist()
            cancel = pass_gzip_format(single_reads)

        # check long reads
        if FQ_HEADERS["LONG"] in samples_df.columns:
            long_reads = samples_df.loc[sample_id, :, :, :][FQ_HEADERS["LONG"]].tolist()
            cancel = pass_gzip_format(long_reads)

        # check short_interleaved
        if FQ_HEADERS["INTERLEAVED"] in samples_df.columns:
            interleaved_reads = samples_df.loc[sample_id, :, :, :][FQ_HEADERS["INTERLEAVED"]].tolist()
            cancel = pass_gzip_format(interleaved_reads)

    if cancel:
        sys.exit(-1)
    else:
        if is_sra:
            return samples_df, "SRA"
        else:
            return samples_df, "FQ"


def parse_mags(mags_dir):
    bin_list = []
    for bin_ in glob.glob(mags_dir + "/*/*.fa"):
        bin_dict = dict()
        bin_dict["path"] = bin_.strip()
        bin_dict["id"] = os.path.basename(bin_).rstrip(".fa")
        bin_list.append(bin_dict)
    mags_df = pd.DataFrame(bin_list).set_index("id", drop=False)
    return mags_df


def get_raw_input_list(wildcards, samples_df, data_type):
    fqs = []
    headers = HEADERS[data_type]
    for k, v in headers.items():
        if v in samples_df.columns:
            reads = get_reads(samples_df, wildcards, v)
            if len(reads) > 0:
                fqs += reads
    return fqs


def get_raw_input_dict(wildcards, samples_df, data_type):
    fqs = {}
    headers = HEADERS[data_type]
    for k, v in headers.items():
        if v in samples_df.columns:
            reads = get_reads(samples_df, wildcards, v)
            if len(reads) > 0:
                fqs[v] = reads
    return fqs


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