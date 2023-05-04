#!/usr/bin/env python

import glob
import os
import sys
import json

from ruamel.yaml import YAML
from executor import execute
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


def get_samples_for_assembly_list(wildcards, samples_df, samples_dir):
    samples_id_list = get_samples_id_by_assembly_and_binning_group(samples_df, wildcards.assembly_group, wildcards.binning_group)
    samples_json_list = []

    for sample_id in samples_id_list:
        jsonfile = os.path.join(samples_dir, f"reads/{sample_id}/{sample_id}.json")
        samples_json_list.append(jsonfile)

    return samples_json_list


def get_samples_for_assembly_dict(wildcards, samples_df, samples_dir):
    samples_dict = {
        "PE_FORWARD": [],
        "PE_REVERSE": [],
        "SE": []
    }

    input_list = get_samples_for_assembly_list(wildcards, samples_df, samples_dir)
    for sample_json in input_list:
        with open(sample_json, "rt") as ih:
            jsondata = json.load(ih)

            if len(jsondata.get("PE_FORWARD", [])) > 0:
                samples_dict["PE_FORWARD"].append(jsondata["PE_FORWARD"])
                samples_dict["PE_REVERSE"].append(jsondata["PE_REVERSE"])

            if len(jsondata.get("SE", [])) > 0:
                samples_dict["SE"].append(jsondata["SE"])

    return samples_dict


def get_samples_for_assembly_megahit(wildcards, samples_df, samples_dir):
    samples_dict = get_samples_for_assembly_dict(wildcards, samples_df, samples_dir)

    reads = ""
    if len(samples_dict.get("PE_FORWARD", [])) > 0:
        reads = f''' -1 {",".join(samples_dict["PE_FORWARD"])} -2 {",".join(samples_dict["PE_REVERSE"])} '''
    if len(samples_dict.get("SE", [])) > 0:
        reads += f''' -r {",".join(samples_dict["SE"])} '''

    return reads


def get_samples_for_assembly_idba_ud(wildcards, samples_df, samples_dir):
    samples_dict = get_samples_for_assembly_dict(wildcards, samples_df, samples_dir)
    reads = os.path.join(
        samples_dir,
        "reads_asm",
        f"{wildcards.binning_group}.{wildcards.assembly_group}/idba_ud/reads_temp")
    reads_dir = os.path.dirname(reads)

    cmd = f'''(mkdir -p {reads_dir})'''

    input_str = ""

    if len(samples_dict.get("PE_FORWARD", [])) > 0:
        pe_forward_list = samples_dict["PE_FORWARD"]
        pe_reverse_list = samples_dict["PE_REVERSE"]
        reads_short = f"{reads}.short.fa"
        if len(pe_forward_list) == 1:
            cmd = f''' && (seqtk mergepe {pe_forward_list[0]} {pe_reverse_list[0]} | seqtk seq -A - > {reads_short})'''
        else:
            cmd = f''' && (cat {" ".join(pe_forward_list)} > {reads}.pe.1.fq.gz)'''
            cmd += f''' && (cat {" ".join(pe_reverse_list)} > {reads}.pe.2.fq.gz)'''
            cmd += f''' && (seqtk mergepe {reads}.pe.1.fq.gz {reads}.pe.2.fq.gz | seqtk seq -A - > {reads_short})'''
            cmd += f''' && (rm -rf {reads}.pe.1.fq.gz {reads}.pe.2.fq.gz)'''
        input_str = f" -r {reads_short} "

    # FIXME can we combined PE and SE data?
    if len(samples_dict.get("SE", [])) > 0:
        se_list = samples_dict["SE"]
        reads_short = f"{reads}.short.fa"
        if len(se_list) == 1:
            cmd += f''' && (seqtk seq -A {se_list[0]} >> {reads_short})'''
        else:
            cmd = f''' && (cat {" ".join(se_list)} > {reads}.se.fq.gz)'''
            cmd += f''' && (seqtk seq -A {reads}.se.fq.gz >> {reads_short})'''
            cmd += f''' && (rm -rf {reads}.se.fq.gz)'''
        input_str = f" -r {reads_short} "

    if len(samples_dict.get("LONG", [])) > 0:
        long_list = samples_dict["LONG"]
        reads_long = f"{reads}.long.fa"
        if len(long_list) == 1:
            cmd += f''' && (seqtk seq -A {long_list[0]} > {reads_long})'''
        else:
            cmd = f''' && (cat {" ".join(long_list)} > {reads}.long.fq.gz)'''
            cmd += f''' && (seqtk seq -A {reads}.long.fq.gz > {reads_long})'''
            cmd += f''' && (rm -rf {reads}.long.fq.gz)'''
        input_str = f" {input_str} -l {reads_long} "

    return cmd, input_str, reads_dir


def get_samples_for_assembly_spades(wildcards, samples_df, samples_dir, assembler):
    """
    ### Read-pair libraries

    By using command line interface, you can specify up to nine different paired-end libraries,
    up to nine mate-pair libraries and also up to nine high-quality mate-pair ones.
    If you wish to use more, you can use [YAML data set file](#yaml).
    We further refer to paired-end and mate-pair libraries simply as to read-pair libraries.

    By default, SPAdes assumes that paired-end and high-quality mate-pair reads have forward-reverse (fr)
    orientation and usual mate-pairs have reverse-forward (rf) orientation.
    However, different orientations can be set for any library by using SPAdes options.

    To distinguish reads in pairs we refer to them as left and right reads.
    For forward-reverse orientation, the forward reads correspond to the left reads and the reverse reads, to the right.
    Similarly, in reverse-forward orientation left and right reads correspond to reverse and forward reads, respectively, etc.


    [
        {
            orientation: "fr",
            type: "paired-end",
            right reads: [
                "/FULL_PATH_TO_DATASET/lib_pe1_right_1.fastq",
                "/FULL_PATH_TO_DATASET/lib_pe1_right_2.fastq"
            ],
            left reads: [
                "/FULL_PATH_TO_DATASET/lib_pe1_left_1.fastq",
                "/FULL_PATH_TO_DATASET/lib_pe1_left_2.fastq"
            ]
        },
        {
            orientation: "rf",
            type: "mate-pairs",
            right reads: [
                "/FULL_PATH_TO_DATASET/lib_mp1_right.fastq"
            ],
            left reads: [
                "/FULL_PATH_TO_DATASET/lib_mp1_left.fastq"
            ]
        },
        {
            type: "single",
            single reads: [
                "/FULL_PATH_TO_DATASET/pacbio_ccs.fastq"
            ]
        },
        {
            type: "pacbio",
            single reads: [
                "/FULL_PATH_TO_DATASET/pacbio_clr.fastq"
            ]
        }
    ]
    """

    yaml = YAML()
    yaml.default_flow_style = False

    samples_dict = get_samples_for_assembly_dict(wildcards, samples_df, samples_dir)
    datasets = []

    have_pe = False
    if len(samples_dict.get("PE_FORWARD", [])) > 0:
        have_pe = True
        datasets_pe = {}
        #datasets_pe["orientation"] = "fr"
        datasets_pe["type"] = "paired-end"
        datasets_pe["left reads"] = [os.path.realpath(fq) for fq in samples_dict["PE_FORWARD"]]
        datasets_pe["right reads"] = [os.path.realpath(fq) for fq in samples_dict["PE_REVERSE"]]
        datasets.append(datasets_pe)

    if assembler == "metaspades":
        if not have_pe:
            print("Currently metaSPAdes supports only a single short-read library which has to be paired-end.")
            sys.exit(-1)

    if assembler != "metaspades":
        if len(samples_dict.get("SE", [])) > 0:
            datasets_se = {}
            datasets_se["type"] = "single"
            datasets_se["single reads"] = [os.path.realpath(fq) for fq in samples_dict["SE"]]
            datasets.append(datasets_se)

    if len(samples_dict.get("LONG", [])) > 0:
        datasets_long = {}
        # FIXME or nanopore ?
        datasets_long["type"] = "pacbio" # or nanopore ? # FIXME
        datasets_long["single reads"] = [os.path.realpath(fq) for fq in samples_dict["LONG"]]

    datasets_yaml = os.path.join(
        samples_dir,
        "reads_asm",
        f"{wildcards.binning_group}.{wildcards.assembly_group}/spades/datasets.yaml")
    os.makedirs(os.path.dirname(datasets_yaml), exist_ok=True)

    with open(datasets_yaml, "wt") as oh:
        yaml.dump(datasets, oh)

    return datasets_yaml


def get_samples_for_assembly_plass(wildcards, samples_df, samples_dir):
    samples_dict = get_samples_for_assembly_dict(wildcards, samples_df, samples_dir)

    reads = os.path.join(
        samples_dir,
        "reads_asm",
        f"{wildcards.binning_group}.{wildcards.assembly_group}/plass/reads_temp")
    reads_dir = os.path.dirname(reads)

    cmd = f'''(mkdir -p {reads_dir})'''

    reads_pe_list = []
    if len(samples_dict.get("PE_FORWARD", [])) > 0:
        for r1, r2 in zip(samples_dict["PE_FORWARD"], samples_dict["PE_REVERSE"]):
            reads_pe_list.append(r1)
            reads_pe_list.append(r2)
    reads_pe = ""
    if len(reads_pe_list) > 0:
        reads_pe = " ".join(reads_pe_list)

    reads_se = ""
    reads_se_list = samples_dict.get("SE", [])
    if len(reads_se_list) == 1:
        reads_se = reads_se_list[0]
    elif len(reads_se_list) > 1:
        reads_se = f"{reads}.se.fq.gz"
        cmd += f''' && (cat {" ".join(reads_se_list)} > {reads_se})'''

    return cmd, reads_pe, reads_se, reads_dir


def get_samples_for_assembly_opera_ms(wildcards, samples_df, samples_dir, asm_dir):
    samples_dict = get_samples_for_assembly_dict(wildcards, samples_df, samples_dir)
    prefix = f"{wildcards.binning_group}.{wildcards.assembly_group}"

    reads = os.path.join(
        samples_dir,
        "reads_asm",
        f"{prefix}/opera_ms/reads_temp")
    reads_dir = os.path.dirname(reads)

    cmd = f'''(mkdir -p {reads_dir})'''

    input_str = ""

    reads_pe1 = ""
    reads_pe2 = ""
    if len(samples_dict.get("PE_FORWARD", [])) > 0:
        pe_forward_list = samples_dict["PE_FORWARD"]
        pe_reverse_list = samples_dict["PE_REVERSE"]
        if len(pe_forward_list) == 1:
            reads_pe1 = pe_forward_list[0]
            reads_pe2 = pe_reverse_list[0]
        else:
            reads_pe1 = f'''{reads}.pe.1.fq.gz'''
            reads_pe2 = f'''{reads}.pe.2.fq.gz'''
            cmd += f''' && (cat {" ".join(pe_forward_list)} > {reads_pe1})'''
            cmd += f''' && (cat {" ".join(pe_reverse_list)} > {reads_pe2})'''
        input_str = f" {input_str} --short-read1 {reads_pe1} --short-read2 {reads_pe2} "

    reads_long = ""
    if len(samples_dict.get("LONG", [])) > 0:
        long_list = samples_dict["LONG"]
        if len(long_list) == 1:
            reads_long = long_list[0]
        else:
            reads_long = f'''{reads}.long.fq.gz'''
            cmd += f''' && (cat {" ".join(long_list)} > {reads_long})'''
        input_str = f" {input_str} --long-read {reads_long} "

    scaftigs_gz = os.path.join(
        asm_dir,
        "scaftigs",
        f"{prefix}.{wildcards.assembler}/{prefix}.{wildcards.assembler}.scaftigs.fa.gz")
    scaftigs = os.path.join(
        samples_dir,
        "reads_asm",
        f"{prefix}/opera_ms/scafitgs_temp.fa")
    cmd += f''' && (pigz -f -dkc {scaftigs_gz} > {scaftigs})'''
    input_str = f" {input_str} --contig-file {scaftigs} "

    return cmd, input_str, reads_dir


def get_samples_for_metaquast(wildcards, samples_df, samples_dir):
    samples_dict = get_samples_for_assembly_dict(wildcards, samples_df, samples_dir)
    prefix = f"{wildcards.binning_group}.{wildcards.assembly_group}.{wildcards.assembler}"

    reads = os.path.join(
        samples_dir,
        "reads_asm",
        f"{prefix}/metaquast/reads_temp")
    reads_dir = os.path.dirname(reads)

    cmd = f'''(mkdir -p {reads_dir})'''

    input_str = ""

    reads_pe1 = ""
    reads_pe2 = ""
    if len(samples_dict.get("PE_FORWARD", [])) > 0:
        pe_forward_list = samples_dict["PE_FORWARD"]
        pe_reverse_list = samples_dict["PE_REVERSE"]
        if len(pe_forward_list) == 1:
            reads_pe1 = pe_forward_list[0]
            reads_pe2 = pe_reverse_list[0]
        else:
            reads_pe1 = f'''{reads}.pe.1.fq.gz'''
            reads_pe2 = f'''{reads}.pe.2.fq.gz'''
            cmd += f''' && (cat {" ".join(pe_forward_list)} > {reads_pe1})'''
            cmd += f''' && (cat {" ".join(pe_reverse_list)} > {reads_pe2})'''
        input_str = f" {input_str} --pe1 {reads_pe1} --pe2 {reads_pe2} "

    reads_se = ""
    if len(samples_dict.get("SE", [])) > 0:
        reads_se_list = samples_dict.get("SE", [])
        if len(reads_se_list) == 1:
            reads_se = reads_se_list[0]
        elif len(reads_se_list) > 1:
            reads_se = f"{reads}.se.fq.gz"
            cmd += f''' && (cat {" ".join(reads_se_list)} > {reads_se})'''
        input_str = f" {input_str} --single {reads_se} "

    reads_long = ""
    if len(samples_dict.get("LONG", [])) > 0:
        long_list = samples_dict["LONG"]
        if len(long_list) == 1:
            reads_long = long_list[0]
        else:
            reads_long = f'''{reads}.long.fq.gz'''
            cmd += f''' && (cat {" ".join(long_list)} > {reads_long})'''
        input_str = f" {input_str} --pacbio {reads_long} "
        # FIXME
        #input_str = f"{input_str} --nanopore {reads_long}"

    return cmd, input_str, reads_dir
