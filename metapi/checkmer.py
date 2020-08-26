#!/usr/bin/env python3

import pandas as pd
import concurrent.futures
import re
import os
import sys
import argparse


def MIMAG_quality_level(row):
    """
    https://doi.org/10.1038/nbt.3893
    """
    if (row["completeness"] > 90.0) and (row["contamination"] < 5.0):
        return "high_quality"
    elif (row["completeness"] > 50.0) and (row["contamination"] < 10.0):
        return "medium_quality"
    else:
        return "low_quality"


def SGB_quality_level(row):
    """
    https://doi.org/10.1016/j.cell.2019.01.001
    """
    if (
        (row["strain_heterogeneity"] < 0.5)
        and (row["completeness"] > 90.0)
        and (row["contamination"] < 5.0)
    ):
        return "high_quality"
    elif (row["completeness"] > 50.0) and (row["contamination"] < 5.0):
        return "medium_quality"
    else:
        return "low_quality"


def quality_score(row):
    """
    https://doi.org/10.1038/s41586-019-0965-1
    """
    return row["completeness"] - 5 * row["contamination"]


def parse(checkm_table):
    if os.path.getsize(checkm_table) > 0:
        checkm_df = pd.read_csv(checkm_table, sep="\t")
        return checkm_df
    else:
        return None


def checkm_report(checkm_list, output, threads):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for df in executor.map(parse, checkm_list):
            if df is not None:
                df_list.append(df)

    df_ = pd.concat(df_list).rename(
        columns={
            "Bin Id": "bin_id",
            "Marker lineage": "marker_lineage",
            "# genomes": "genomes",
            "# markers": "markers",
            "# marker sets": "marker_sets",
            "Completeness": "completeness",
            "Contamination": "contamination",
            "Strain heterogeneity": "strain_heterogeneity",
        }
    )

    df_["MIMAG_quality_level"] = df_.apply(lambda x: MIMAG_quality_level(x), axis=1)
    df_["SGB_quality_level"] = df_.apply(lambda x: SGB_quality_level(x), axis=1)
    df_["quality_score"] = df_.apply(lambda x: quality_score(x), axis=1)

    df_.to_csv(output, sep="\t", index=False)
    return df_


def main():
    parser = argparse.ArgumentParser("CheckM reporter")
    parser.add_argument("--checkm_list", type=str, help="checkm out list")
    parser.add_argument("--output", type=str, required=True, help="checkm output file")
    parser.add_argument(
        "--threads", type=int, default=8, help="threads used on combine CheckM output"
    )
    args = parser.parse_args()

    checkm_list = [l.strip() for l in open(args.checkm_list, "r")]
    checkm_report(checkm_list, args.output, args.threads)


if __name__ == "__main__":
    main()
