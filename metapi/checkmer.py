#!/usr/bin/env python3

import pandas as pd
import concurrent.futures
import re
import os
import sys
import argparse


def MIMAG_quality_level(row):
    """
    https://www.nature.com/articles/nbt.3893
    """
    if (row["completeness"] > 90.0) and (row["contamination"] < 5.0):
        return "high_quality"
    elif (row["completeness"] > 50.0) and (row["contamination"] < 10.0):
        return "medium_quality"
    else:
        return "low_quality"


def SGB_quality_level(row):
    """
    http://www.nature.com/articles/s41586-019-0965-1
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


def parse(checkm_file):
    total_list = []
    with open(checkm_file, "r") as ih:
        next(ih), next(ih), next(ih)
        for line in ih:
            if not line.startswith("--"):
                x_x = re.split(r"\s+", line.strip())
                total_list.append(
                    [x_x[0]]
                    + ["-".join(x_x[1:3])]
                    + [
                        int(x_x[3]),
                        int(x_x[4]),
                        int(x_x[5]),
                        int(x_x[6]),
                        int(x_x[7]),
                        int(x_x[8]),
                        int(x_x[9]),
                        int(x_x[10]),
                        int(x_x[11]),
                    ]
                    + [float(x_x[12]), float(x_x[13]), float(x_x[14])]
                )

    checkm_df = pd.DataFrame(
        total_list,
        columns=[
            "bin_id",
            "marker_lineage",
            "genomes",
            "markers",
            "marker_sets",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5+",
            "completeness",
            "contamination",
            "strain_heterogeneity",
        ],
    )
    return checkm_df


def report(checkm_list, output, threads):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_works=threads) as executor:
        for df in executor.map(parse, checkm_list):
            if not df.empty():
                df_list.append(df)

    df_ = pd.concat(df_list)

    df_["MIMAG_quality_level"] = df_.apply(lambda x: MIMAG_quality_level(x), axis=1)
    df_["SGB_quality_level"] = df_.apply(lambda x: SGB_quality_level(x), axis=1)

    df_.to_csv(output, sep='\t')


def main():
    parser = argparse.ArgumentParser("CheckM reporter")
    parser.add_argument("--checkm_list", type=str, help="checkm out list")
    parser.add_argument(
        "--output", type=str, required=True, help="checkm output file"
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="threads used on combine CheckM output"
    )
    args = parser.parse_args()

    checkm_list = [l.strip() for l in open(args.checkm_list, 'r')]
    report(checkm_list, args.output, args.threads)


if __name__ == '__main__':
    main()
