#!/usr/bin/env python

import pandas as pd
import concurrent.futures
import os
import sys
import argparse

def global_init(index_metadata):
    global INDEX_METADATA__
    INDEX_METADATA__ = pd.read_csv(index_metadata, sep='\t')


def get_mgs_id(row):
    return "_".join(row["ID"].split("_")[0:-1])


def get_abun_df_hsx(abun_file):
    sample_id = os.path.basename(abun_file).split(".")[0]

    try:
        if os.path.exists(abun_file):
            abun = pd.read_csv(abun_file, sep='\t')
        else:
            print("%s is not exists" % abun_file)
            return None, None
    except pd.io.common.EmptyDataError:
        print("%s is empty" % abun_file)
        return None, None

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


def get_abun_df_jgi(depth_file):
    sample_id = os.path.basename(depth_file).split(".")[0]

    try:
        if os.path.exists(depth_file):
            depth = pd.read_csv(depth_file, sep='\t')
        else:
            print("%s is not exists" % depth_file)
            return None, None
    except pd.io.common.EmptyDataError:
        print("%s is empty" % depth_file)
        return None, None

    depth = depth.rename(columns={"contigName": "contig_name"})\
                 .merge(INDEX_METADATA__)\
                 .groupby("mgs_id")\
                 .agg({"totalAvgDepth": "mean"})
    depth[sample_id] = depth["totalAvgDepth"] / sum(depth["totalAvgDepth"])
    depth_df = depth.loc[:, ["totalAvgDepth"]].rename(columns={"totalAvgDepth": sample_id})
    abun_df = depth.loc[:, [sample_id]]
    return depth_df, abun_df


def get_all_abun_df(abun_files, workers, func):
    count_list = []
    abun_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for count_df, abun_df in executor.map(func, abun_files):
            if (not count_df is None) and (not abun_df is None):
                count_list.append(count_df)
                abun_list.append(abun_df)

    count_df_ = pd.concat(count_list, axis=1).reset_index()
    abun_df_ = pd.concat(abun_list, axis=1).reset_index()

    return count_df_, abun_df_


def main():
    parser = argparse.ArgumentParser('merge many samples abundance file to one profile')
    parser.add_argument(
        '-l',
        '--abundance_list',
        type=str,
        help='abundance list')
    parser.add_argument(
        '--method',
        default="hsx",
        choices=["hsx", "jgi"],
        help='compute method'
    )
    parser.add_argument(
        '--database',
        default=None,
        help='contig and genome relationships'
    )
    parser.add_argument(
        '--threads',
        default=8,
        type=int,
        help='threads'
    )
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

    if args.method == "jgi" and args.database is None:
        print("pleas supply database when parse jgi depth file")
        sys.exit(1)

    if args.method == "hsx":
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, get_abun_df_hsx)
    elif args.method == "jgi":
        global_init(args.database)
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, get_abun_df_jgi)
    else:
        print("unsupport method: %s" % args.method)

    count_df.to_csv(args.out_count_profile, sep='\t', index=False)
    abun_df.to_csv(args.out_abundance_profile, sep='\t', index=False)


if __name__ == '__main__':
    main()
