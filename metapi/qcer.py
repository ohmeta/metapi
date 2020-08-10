#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
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

    if "output" in kwargs:
        df.to_csv(kwargs["output"], sep="\t", index=False)
    return df


def qc_bar_plot(df, engine, stacked=False, **kwargs):
    if engine == "seaborn":
        # seaborn don't like stacked barplot
        f, ax = plt.subplots(figsize=(10, 7))
        df_ = df.query('reads=="fq1"')
        sns.barplot(x="id", y="num_seqs", hue="step", data=df_)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=-90)

    elif engine == "pandas":
        if not stacked:
            df_ = (
                df.query('reads=="fq1"')
                .pivot(index="id", columns="step", values="num_seqs")
                .loc[:, ["raw", "trimming", "rmhost"]]
            )
            df_.plot(kind="bar", figsize=(10, 7))

        else:
            dict_ = {"id": [], "clean": [], "rmhost": [], "trim": []}
            df = df.set_index("id")
            for i in df.index.unique():
                reads_total = 0

                reads_trimmed = 0
                reads_host = 0
                reads_clean = 0

                if not df.loc[i,].query('reads=="fq1" and step=="raw"').empty:
                    reads_total = df.loc[i,].query('reads=="fq1" and step=="raw"')[
                        "num_seqs"
                    ][0]

                if not df.loc[i,].query('reads=="fq1" and step=="trimming"').empty:
                    reads_trim = df.loc[i,].query('reads=="fq1" and step=="trimming"')[
                        "num_seqs"
                    ][0]

                if not df.loc[i,].query('reads=="fq1" and step=="rmhost"').empty:
                    reads_clean = df.loc[i,].query('reads=="fq1" and step=="rmhost"')[
                        "num_seqs"
                    ][0]

                reads_trimmed = reads_total - reads_trim
                reads_host = reads_trim - reads_clean

                dict_["id"].append(i)
                dict_["trim"].append(reads_trimmed)
                dict_["rmhost"].append(reads_host)
                dict_["clean"].append(reads_clean)

            df_ = pd.DataFrame(dict_).sort_values("id").set_index("id")

            colors = ["#2ca02c", "#ff7f0e", "#1f77b4"]

            df_.plot(kind="bar", stacked=True, color=colors, figsize=(10, 7))

    plt.xlabel("Sample ID")
    plt.ylabel("The number of reads(-pair)")
    plt.title("Fastq quality control barplot", fontsize=11)

    if "output" in kwargs:
        plt.savefig(kwargs["output"])


def main():
    parser = argparse.ArgumentParser(description="quality control reporter")
    parser.add_argument("--raw_stats_list", help="raw stats list")
    parser.add_argument("--trimming_stats_list", help="trimming stats list")
    parser.add_argument("--rmhost_stats_list", help="rmhost stats list")
    parser.add_argument("--output", help="quality control output basename")
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
        tooler.parse,
        8,
    )

    df_ = compute_host_rate(df, output=os.path.join(args.output, ".stats.tsv"))
    qc_bar_plot(df_, "seaborn", output=os.path.join(args.output, ".plot.pdf"))


if __name__ == "__main__":
    main()
