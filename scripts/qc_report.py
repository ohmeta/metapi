#!/usr/bin/env python3

import pandas as pd
import concurrent.futures
import subprocess
import argparse
import os

def get_reads(df, id_, col_):
    return df.loc[[id_], col_].dropna().tolist()


def run_(tuple_):
    try:
        output = subprocess.check_output(
            tuple_[0], shell=True, stderr=subprocess.STDOUT, universal_newlines=True
        )
    except subprocess.CalledProcessError as e:
        print(e.output)
        return None

    out_list = output.split("\n")
    header = out_list[0].split("\t")
    data = []

    for line in out_list[1:]:
        content = tuple(line.split("\t"))
        data.append(content)

    df = pd.DataFrame(data, columns=header)
    df["id"] = tuple_[1]
    df["step"] = tuple_[2]
    df["fq_type"] = tuple_[3]
    df["reads"] = tuple_[4]
    return df


def run(cmd_list, workers):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(run_, cmd_list):
            if df is not None:
                df_list.append(df)
    df_ = pd.concat(df_list)
    return df_


def gen(fastq_list, step, fq_encoding, is_pe=True):
    fq_df = pd.read_csv(fastq_list, sep="\t").set_index("id")
    cmd_list = []

    for i in fq_df.index.unique():
        fq1_list = get_reads(fq_df, i, "fq1")
        if is_pe:
            fq2_list = get_reads(fq_df, i, "fq2")

        if is_pe:
            if len(fq1_list) == 1:
                cmd = "seqkit stats -a -T -b -j 1 -E %s %s %s" % (
                    fq_encoding,
                    fq1_list[0],
                    fq2_list[0],
                )
                cmd_list.append((cmd, i, step, "pe", ["fq1", "fq2"]))
            else:
                cmd_1 = "cat %s | seqkit stats -a -T -b -j 1 -E %s" % (
                    " ".join(fq1_list),
                    fq_encoding,
                )
                cmd_list.append((cmd_1, i, step, "pe", ["fq1"]))
                cmd_2 = "cat %s | seqkit stats -a -T -b -j 1 -E %s" % (
                    " ".join(fq2_list),
                    fq_encoding,
                )
                cmd_list.append((cmd_2, i, step, "pe", ["fq2"]))
        else:
            cmd = "cat %s | seqkit stats -a -T -b -j 1 -E %s" % (
                " ".join(fq1_list),
                fq_encoding,
            )
            cmd_list.append((cmd, i, step, "se", ["fq1"]))
    return cmd_list


def main():
    parser = argparse.ArgumentParser(description="generate quality control report from raw, trimming, rmhost data")
    parser.add_argument("--raw_list", help="raw data list, headers: id fq1 fq2")
    parser.add_argument("--trimming_list", help="trimming data list, headers: id fq1 fq2")
    parser.add_argument("--rmhost_list", help="rmhost data list, headers: id fq1 fq2")
    parser.add_argument("--is_se", action='store_true', default=False, help='default: is_pe')
    parser.add_argument("--fq_encoding", help="fastq quality encoding, default: sanger", default="sanger")
    parser.add_argument("--threads", help="threads, default: 8", default=8)
    parser.add_argument("--output", help="qc report output")

    args = parser.parse_args()

    cmd_raw = gen(args.raw_list, "raw", args.fq_encoding, not args.is_se)
    cmd_trimming = gen(args.trimming_list, "trimming", args.fq_encoding, not args.is_se)
    cmd_rmhost = gen(args.rmhost_list, "rmhost", args.fq_encoding, not args.is_se)

    cmd = cmd_raw + cmd_trimming + cmd_rmhost

    df = run(cmd, args.threads)
    df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
