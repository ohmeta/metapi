#!/usr/bin/env python

import pandas as pd
import os
from concurrent import futures
import argparse

def get_mgs_id(row):
    return "_".join(row["ID"].split("_")[0:-1])


def get_abun_df(abun_file):
    sample_id = os.path.basename(abun_file).split(".")[0]
    abun = pd.read_csv(abun_file, sep='\t')
    abun["mgs_id"] = abun.apply(get_mgs_id, axis=1)

    count_df = abun.loc[:, ["mgs_id", "reads_pairs"]]\
                   .groupby("mgs_id")\
                   .agg({"reads_pairs": 'sum'})\
                   .rename(columns={"reads_pairs": sample_id})
    abun_df = abun.loc[:, ["mgs_id", "gene_abundance"]]\
                  .groupby("mgs_id")\
                  .agg({"gene_abundance": 'sum'})\
                  .rename(columns={"gene_abundance": sample_id})
    return count_df, abun_df


def get_all_abun_df(abun_files):
    count_list = []
    abun_list = []
    with futures.ProcessPoolExecutor() as pool:
        for count_df, abun_df in pool.map(get_abun_df, abun_files):
            count_list.append(count_df)
            abun_list.append(abun_df)

    count_df_ = pd.concat(count_list, axis=1)
    abun_df_ = pd.concat(abun_list, axis=1)

    return count_df_, abun_df_


def main():
    parser = argparse.ArgumentParser('merge many samples abundance file to one profile')
    parser.add_argument(
        '-l',
        '--abundance_list',
        type=str,
        help='abundance list')
    parser.add_argument(
        '--out_count_profile',
        type=str,
        help='output count profile')
    parser.add_argument(
        '--out_abundance_profile',
        type=str,
        help='output abundance profile')
    args = parser.parse_args()

    abun_files = pd.read_csv(args.abundance_list, names=["path"])\
                   .loc[:, "path"].values

    count_df, abun_df = get_all_abun_df(abun_files)

    count_df.reset_index().to_csv(args.out_count_profile, sep='\t', index=False)
    abun_df.reset_index().to_csv(args.out_abun_profile, sep='\t', index=False)
