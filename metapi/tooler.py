#!/usr/bin/env python3

import os
import concurrent.futures
import pandas as pd


def parse(stats_file):
    try:
        if os.path.exists(stats_file):
            df = pd.read_csv(stats_file, sep="\t")
            if not df.empty:
                return df
            else:
                return None
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

    if "output" in kwargs:
        df_.to_csv(kwargs["output"], sep="\t", index=False)
    return df_


def merge2(input_list, func, workers, **kwargs):
    df1_list = []
    df2_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df1, df2 in executor.map(func, input_list):
            if df1 is not None:
                df1_list.append(df1)
            if df2 is not None:
                df2_list.append(df2)

    df_1 = pd.concat(df1_list)
    df_2 = pd.concat(df2_list)

    if "output_1" in kwargs:
        df_1.to_csv(kwargs["output_1"], sep="\t", index=False)
    if "output_2" in kwargs:
        df_2.to_csv(kwargs["output_2"], sep="\t", index=False)

    return df_1, df_2
