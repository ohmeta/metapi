#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from metapi import tooler


def change(output, sample_id, step, fq_type, reads_list):
    df = pd.read_csv(output, sep="\t").sort_values("file", ascending=True)
    df["id"] = sample_id
    df["reads"] = reads_list
    df["step"] = step
    df["fq_type"] = fq_type
    df.to_csv(output, sep="\t", index=False)


def compute_host_rate(df, **kwargs):
    host_rate = {}
    df = df.set_index("id")
    for i in df.index.unique():
        if not df.loc[i,].query('reads=="fq1" and step=="rmhost"').empty:
            reads_num_rmhost = df.loc[i,].query('reads=="fq1" and step=="rmhost"')[
                "num_seqs"
            ][0]
            if not df.loc[i,].query('reads=="fq1" and step=="trimming"').empty:
                reads_num = df.loc[i,].query('reads=="fq1" and step=="trimming"')[
                    "num_seqs"
                ][0]
            elif not df.loc[i,].query('reads=="fq1" and step=="raw"').empty:
                reads_num = df.loc[i,].query('reads=="fq1" and step=="raw"')[
                    "num_seqs"
                ][0]
            hostrate = (reads_num - reads_num_rmhost) / reads_num
            host_rate[i] = hostrate
        else:
            host_rate[i] = np.nan

    df = df.reset_index()
    df["host_rate"] = df.apply(lambda x: host_rate[x["id"]], axis=1)

    if "save" in kwargs:
        if kwargs["save"]:
            if "output" in kwargs:
                df.to_csv(kwargs["output"], sep="\t", index=False)
            else:
                print("please specific output parameter")
    return df


def main():
    parser = argparse.ArgumentParser(description="quality control reporter")
    parser.add_argument("--raw_stats_list", help="raw stats list")
    parser.add_argument("--trimming_stats_list", help="trimming stats list")
    parser.add_argument("--rmhost_stats_list", help="rmhost stats list")
    parser.add_argument("--output", help="quality control output")
    args = parser.parse_args()

    raw_list = pd.read_csv(args.raw_stats_list, header=None, names=["raw"])
    trimming_list = pd.read_csv(
        args.trimming_stats_list, header=None, names=["trimming"]
    )
    rmhost_list = pd.read_csv(args.rmhost_stats_list, header=None, names=["rmhost"])

    df = tooler.merge(
        raw_list["raw"].dropna().tolist()
        + trimming_list["trimming"].dropna().tolist()
        + rmhost_list["rmhost"].dropna().tolist(),
        8,
    )

    compute_host_rate(df, save=True, output=args.output)


if __name__ == "__main__":
    main()
