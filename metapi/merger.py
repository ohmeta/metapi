#!/usr/bin/env python

import pandas as pd
import concurrent.futures
import os
import sys
import argparse


def global_init(index_metadata):
    global INDEX_METADATA__
    INDEX_METADATA__ = pd.read_csv(index_metadata, sep='\t')


def set_lineages_to(row, key, level):
    LINEAGES = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
    LINEAGES_DICT = {
        "strain": LINEAGES,
        "species": LINEAGES[0:7],
        "genus": LINEAGES[0:6],
        "family": LINEAGES[0:5],
        "order": LINEAGES[0:4],
        "class": LINEAGES[0:3],
        "phylum": LINEAGES[0:2],
        "superkingdom": LINEAGES[0:1],
    }

    LEVEL_DICT = {
        "strain": "t",
        "species": "s",
        "genus": "g",
        "family": "f",
        "order": "o",
        "class": "c",
        "phylum": "p",
        "superkingdom": "k"
    }

    lineages_dict = {
        "k": "k__unclassified_" + row["mgs_id"],
        "p": "p__unclassified_" + row["mgs_id"],
        "c": "c__unclassified_" + row["mgs_id"],
        "o": "o__unclassified_" + row["mgs_id"],
        "f": "f__unclassified_" + row["mgs_id"],
        "g": "g__unclassified_" + row["mgs_id"],
        "s": "s__unclassified_" + row["mgs_id"],
        "t": "t__unclassified_" + row["mgs_id"],
    }
    for line in row[key].split(";"):
        lev, tax = line.split("__")
        if lev == "d":
            lev = "k"
        if tax != "":
            lineages_dict[lev] = lev + "__" + tax
            if lev == "s" or lev == "t":
                lineages_dict[lev] = lineages_dict[lev] + "_" + row["mgs_id"]

    lineages = []
    for i in LINEAGES_DICT[level]:
        lineages.append(lineages_dict[LEVEL_DICT[i]])

    return "|".join(lineages)


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


def get_all_abun_df(abun_files, workers, method):
    count_list = []
    abun_list = []
    if method == "jgi":
        func = get_abun_df_jgi
    elif method == "hsx":
        func = get_abun_df_hsx
    else:
        print("unspoort method %s" % method)

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for count_df, abun_df in executor.map(func, abun_files):
            if (not count_df is None) and (not abun_df is None):
                count_list.append(count_df)
                abun_list.append(abun_df)

    count_df_ = pd.concat(count_list, axis=1).reset_index()
    abun_df_ = pd.concat(abun_list, axis=1).reset_index()

    return count_df_, abun_df_


def get_profile(abun_tax_df, samples_list, key, profile_tsv):
    #level_ = "lineages_" + level + "_new"
    #abun_tax_df[level_] = abun_tax_df.apply(lambda x: set_lineages_to(x, level), axis=1)
    level_ = key
    _profile = abun_tax_df.loc[:, [level_] + samples_list]\
                       .melt(id_vars=[level_])\
                       .groupby(["variable", level_])\
                       .agg({"value": 'sum'})\
                       .reset_index()\
                       .pivot(index=level_,
                              columns="variable",
                              values="value")
    profile_ = _profile.reset_index().loc[:, [level_] + samples_list]
    profile_.to_csv(profile_tsv, sep='\t', index=False)


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
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, "hsx")
    elif args.method == "jgi":
        global_init(args.database)
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, "jgi")
    else:
        print("unsupport method: %s" % args.method)

    count_df.to_csv(args.out_count_profile, sep='\t', index=False)
    abun_df.to_csv(args.out_abundance_profile, sep='\t', index=False)


if __name__ == '__main__':
    main()
