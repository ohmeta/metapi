#!/usr/bin/env python3

import pandas as pd
import concurrent.futures
import os


def change(output, sample_id, step, fq_type, reads_list):
    df = pd.read_csv(output, sep='\t').sort_values("file", ascending=True)
    df["id"] = sample_id
    df["reads"] = reads_list
    df["step"] = step
    df["fq_type"] = fq_type
    df.to_csv(output, sep='\t', index=False)


def parse(stats_file):
    try:
        if os.path.exists(stats_file):
            df = pd.read_csv(stats_file, sep='\t')
            return df
        else:
            print("%s is not exists" % stats_file)
            return None
    except pd.io.common.EmptyDataError:
        print("%s is empty, please check" % stats_file)
        return None


def merge(input_list, output, workers):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(parse, input_list):
            if df is not None:
                df_list.append(df)

    df_ = pd.concat(df_list)
    df_.to_csv(output, sep='\t', index=False)
