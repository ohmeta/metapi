#!/usr/bin/env python3

import os
import concurrent.futures
import pandas as pd


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
