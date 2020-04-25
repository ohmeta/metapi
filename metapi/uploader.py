#!/usr/bin/env python3

import os
import pandas as pd
import concurrent.futures


def gen_samples_info(samples, output, config):
    samples_df = pd.DataFrame(
        columns=["project_accession", "sample_name"]
        + list(config["upload"]["samples"].keys())
    )
    samples_df["sample_name"] = samples.index.unique()
    samples_df["project_accession"] = config["upload"]["project_accession"]
    for key in config["upload"]["samples"]:
        samples_df[key] = config["upload"]["samples"][key]
    samples_df.to_excel(output, index=False)


def parse_md5(md5_file):
    try:
        if os.path.exists(md5_file):
            df = pd.read_csv(
                md5_file, sep="\s+", header=None, names=["file_md5", "file_name"]
            )
            if df.empty:
                print("%s is empty, please check" % md5_file)
                return None
            df["sample_name"] = df.apply(
                lambda x: os.path.basename(x["file_name"]).split(".")[0], axis=1
            )
            df["file_name"] = df.apply(
                lambda x: os.path.basename(x["file_name"]), axis=1
            )
            if len(df) == 2:
                df_fq1 = df.iloc[0].to_frame().T
                df_fq2 = (
                    df.iloc[1]
                    .to_frame()
                    .T.rename(
                        columns={"file_name": "file2_name", "file_md5": "file2_md5"}
                    )
                )
                return pd.merge(df_fq1, df_fq2)
            else:
                df["file2_name"] = ""
                df["file2_md5"] = ""
                return df
        else:
            print("%s is not exists" % md5_file)
            return None
    except pd.io.common.EmptyDataError:
        print("%s is empty, please check" % md5_file)
        return None


def gen_info(input_list, output, config, workers, group):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df_ in executor.map(parse_md5, input_list):
            if df_ is not None:
                df_list.append(df_)
    df = pd.concat(df_list)

    for key in config["upload"][group].keys():
        df[key] = config["upload"][group][key]

    df["project_accession"] = config["upload"]["project_accession"]

    if group == "sequencing_run":
        run_df = df.loc[
            :,
            ["project_accession", "sample_name"]
            + list(config["upload"][group].keys())
            + ["file_name", "file_md5", "file2_name", "file2_md5"],
        ]
        run_df.to_excel(output, sheet_name="Metadata", index=False)

    if group == "assembly":
        asm_df = df.rename(
            columns={"file_name": "fasta_file_name", "file_md5": "fasta_file_md5"}
        )
        asm_df["assembly_name"] = df["sample_name"] + "-" + df["assembly_method"]
        asm_df["sample_accession"] = ""
        asm_df = asm_df.loc[
            :,
            ["project_accession", "sample_accession", "sample_name", "assembly_name"]
            + list(config["upload"][group].keys())
            + ["fasta_file_name", "fasta_file_md5"],
        ]
        asm_df.to_excel(output, sheet_name="Genome_Assembly", index=False)
