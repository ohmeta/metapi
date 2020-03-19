#!/usr/bin/env python3

import pandas as pd
import concurrent.futures
import re
import os
import sys
import argparse


def parse(checkm_file):
    total_list = []
    with open(checkm_file, "r") as ih:
        next(ih), next(ih), next(ih)
        for line in ih:
            if not line.startswith("--"):
                line_list = re.split(r"\s+", line.strip())
                total_list.append(
                    [line_list[0]] + ["-".join(line_list[1:3])] + line_list[3:15]
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


def report(checkm_list, output, completeness, contamination, threads):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_works=threads) as executor:
        for df in executor.map(parse, checkm_list):
            if not df.empty():
                df_list.append(df)

    df_ = pd.concat(df_list)
    df_ = df_.sort_values(
        by=["completeness", "contamination", "strain_heterogeneity"],
        ascending=[False, True, True],
    )
    return df_
