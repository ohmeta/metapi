#!/usr/bin/env python3

import pandas as pd
import numpy as np
import concurrent.futures
import os
import argparse


def change(output, sample_id, step, fq_type, reads_list):
    df = pd.read_csv(output, sep="\t").sort_values("file", ascending=True)
    df["id"] = sample_id
    df["reads"] = reads_list
    df["step"] = step
    df["fq_type"] = fq_type
    df.to_csv(output, sep="\t", index=False)


def parse(stats_file):
    try:
        if os.path.exists(stats_file):
            df = pd.read_csv(stats_file, sep="\t")
            return df
        else:
            print("%s is not exists" % stats_file)
            return None
    except pd.io.common.EmptyDataError:
        print("%s is empty, please check" % stats_file)
        return None


def merge(input_list, func, workers, **kwargs):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(func, input_list):
            if df is not None:
                df_list.append(df)

    df_ = pd.concat(df_list)
    if "save" in kwargs:
        if kwargs["save"]:
            if "output" in kwargs:
                df_.to_csv(kwargs["output"], sep="\t", index=False)
            else:
                print("please specific output parameter")
    return df_


def N50_N80_N90(x):
    len_total = sum(x)
    len_50, len_80, len_90 = len_total * 0.5, len_total * 0.8, len_total * 0.9
    N50, N80, N90, len_sum = 0, 0, 0, 0
    for i in sorted(x)[::-1]:
        len_sum += i
        if (N50 == 0) and (len_sum >= len_50):
            N50 = i
        if (N80 == 0) and (len_sum >= len_80):
            N80 = i
        if (N90 == 0) and (len_sum >= len_90):
            N90 = i
    return N50, N80, N90


def GC_content(x):
    return (x[("#G", "sum")] + x[("#C", "sum")]) / x[("length", "sum")]


def global_init(contigs_length_range):
    global CONTIGS_LENGTH_RANGES__
    temp = sorted(list(set(contigs_length_range)))
    if (temp is None) or (temp == [0]) or (len(temp) == 0):
        CONTIGS_LENGTH_RANGES__ = [(0, 1500), (1500, 2000), (2000, 2500), 2500]
    elif temp[0] == 0:
        CONTIGS_LENGTH_RANGES__ = [(i, j) for i, j in zip(temp[0:-1], temp[1:])] + [
            temp[-1]
        ]
    else:
        CONTIGS_LENGTH_RANGES__ = [(i, j) for i, j in zip([0] + temp[0:-1], temp)] + [
            temp[-1]
        ]


def cal_len_range(x):
    len_dict = {i: 0 for i in CONTIGS_LENGTH_RANGES__}
    for seq_len in x:
        if seq_len >= CONTIGS_LENGTH_RANGES__[-1]:
            len_dict[CONTIGS_LENGTH_RANGES__[-1]] += 1
        else:
            for i in CONTIGS_LENGTH_RANGES__[:-1]:
                if i[0] <= seq_len < i[1]:
                    len_dict[i] += 1
                    break
    return tuple([len_dict[i] for i in CONTIGS_LENGTH_RANGES__])


def parse_assembly(stats_file):
    df_ = parse(stats_file)
    if df_ is not None:
        df = (
            df_.groupby(["sample_id", "assembler"])
            .agg(
                {
                    "chr": ["count"],
                    "length": ["sum", "min", "max", "std", N50_N80_N90, cal_len_range],
                    "#A": ["sum", "min", "max", "std"],
                    "#C": ["sum", "min", "max", "std"],
                    "#G": ["sum", "min", "max", "std"],
                    "#T": ["sum", "min", "max", "std"],
                    "#CpG": ["sum", "min", "max", "std"],
                }
            )
            .reset_index()
        )
        df["GC_content"] = df.apply(lambda x: GC_content(x), axis=1)

        N_ = ["N50", "N80", "N90"]
        for i in range(0, len(N_)):
            df[("length", N_[i])] = df.apply(
                lambda x: x[("length", "N50_N80_N90")][i], axis=1
            )

        for i in range(0, len(CONTIGS_LENGTH_RANGES__) - 1):
            len_tuple = CONTIGS_LENGTH_RANGES__[i]
            len_range_str = "[%d, %d)" % (len_tuple[0], len_tuple[1])
            df[("length", len_range_str)] = df.apply(
                lambda x: x[("length", "cal_len_range")][i], axis=1
            )

        df[("length", "[%d, )" % CONTIGS_LENGTH_RANGES__[-1])] = df.apply(
            lambda x: x[("length", "cal_len_range")][-1], axis=1
        )
        return df
    else:
        return None


def compute_host_rate(df, **kwargs):
    host_rate = {}
    df = df.set_index("id")
    for i in df.index.unique():
        if (not df.loc[i,].query('reads=="fq1" and step=="rmhost"').empty) and (
            not df.loc[i,].query('reads=="fq1" and step=="trimming"').empty
        ):
            reads_num_rmhost = df.loc[i,].query('reads=="fq1" and step=="rmhost"')[
                "num_seqs"
            ][0]
            reads_num_trimming = df.loc[i,].query('reads=="fq1" and step=="trimming"')[
                "num_seqs"
            ][0]

            hostrate = (reads_num_trimming - reads_num_rmhost) / reads_num_trimming
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

    df = merge(
        raw_list["raw"].dropna().tolist()
        + trimming_list["trimming"].dropna().tolist()
        + rmhost_list["rmhost"].dropna().tolist(),
        8,
    )

    compute_host_rate(df, save=True, output=args.output)


if __name__ == "__main__":
    main()
