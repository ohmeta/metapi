#!/usr/bin/env python

import pandas as pd
import concurrent.futures
import os
import sys
import argparse


def metaphlan_init(version):
    global METAPHLAN_VERSION
    METAPHLAN_VERSION = version


def read_metaphlan_table(table):
    sample_id = os.path.basename(table).split(".")[0]

    if METAPHLAN_VERSION == 2:
        dict_ = {"clade_name": [], sample_id: []}
        with open(table, "r") as ih:
            for line in ih:
                if not line.startswith("#"):
                    clade_name, abun = line.strip().split("\t")[0:2]
                    dict_["clade_name"].append(clade_name)
                    dict_[sample_id].append(abun)
        df = pd.DataFrame(dict_).set_index("clade_name")
        return df

    elif METAPHLAN_VERSION == 3:
        dict_ = {"clade_name": [], "clade_taxid": [], sample_id: []}
        with open(table, "r") as ih:
            for line in ih:
                if not line.startswith("#"):
                    clade_name, clade_taxid, abun = line.strip().split("\t")[0:3]
                    dict_["clade_name"].append(clade_name)
                    dict_["clade_taxid"].append(clade_taxid)
                    dict_[sample_id].append(abun)
        df = pd.DataFrame(dict_).set_index(["clade_name", "clade_taxid"])
        return df
    else:
        return None


def merge_metaphlan_tables(table_files, workers, **kwargs):
    abun_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for abun_df in executor.map(read_metaphlan_table, table_files):
            if abun_df is not None:
                abun_list.append(abun_df)

    abun_df_ = pd.concat(abun_list, axis=1).fillna(0)
    if "output" in kwargs:
        abun_df_.reset_index().to_csv(kwargs["output"], sep="\t", index=False)
    return abun_df_


def profiler_init(index_metadata):
    global INDEX_METADATA__
    INDEX_METADATA__ = pd.read_csv(index_metadata, sep="\t")


def set_lineages_to(row, key, level):
    LINEAGES = [
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
    ]
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
        "superkingdom": "k",
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
            abun = pd.read_csv(abun_file, sep="\t")
        else:
            print("%s is not exists" % abun_file)
            return None, None
    except pd.io.common.EmptyDataError:
        print("%s is empty" % abun_file)
        return None, None

    abun["mgs_id"] = abun.apply(get_mgs_id, axis=1)

    count_df = (
        abun.loc[:, ["mgs_id", "reads_pairs"]]
        .groupby("mgs_id")
        .agg({"reads_pairs": "sum"})
        .rename(columns={"reads_pairs": sample_id})
    )
    abun_df = (
        abun.loc[:, ["mgs_id", "gene_abundance"]]
        .groupby("mgs_id")
        .agg({"gene_abundance": "sum"})
        .rename(columns={"gene_abundance": sample_id})
    )
    return count_df, abun_df


def get_abun_df_jgi(depth_file):
    sample_id = os.path.basename(depth_file).split(".")[0]

    try:
        if os.path.exists(depth_file):
            depth = pd.read_csv(depth_file, sep="\t")
        else:
            print("%s is not exists" % depth_file)
            return None, None
    except pd.io.common.EmptyDataError:
        print("%s is empty" % depth_file)
        return None, None

    depth = (
        depth.rename(columns={"contigName": "contig_name"})
        .merge(INDEX_METADATA__)
        .groupby("mgs_id")
        .agg({"totalAvgDepth": "mean"})
    )
    depth[sample_id] = depth["totalAvgDepth"] / sum(depth["totalAvgDepth"])
    depth_df = depth.loc[:, ["totalAvgDepth"]].rename(
        columns={"totalAvgDepth": sample_id}
    )
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
            if (count_df is not None) and (abun_df is not None):
                count_list.append(count_df)
                abun_list.append(abun_df)

    count_df_ = pd.concat(count_list, axis=1).reset_index()
    abun_df_ = pd.concat(abun_list, axis=1).reset_index()

    return count_df_, abun_df_


def get_profile(abun_tax_df, samples_list, key, profile_tsv):
    # level_ = "lineages_" + level + "_new"
    # abun_tax_df[level_] = abun_tax_df.apply(lambda x: set_lineages_to(x, level), axis=1)
    level_ = key
    _profile = (
        abun_tax_df.loc[:, [level_] + samples_list]
        .melt(id_vars=[level_])
        .groupby(["variable", level_])
        .agg({"value": "sum"})
        .reset_index()
        .pivot(index=level_, columns="variable", values="value")
    )
    profile_ = _profile.reset_index().loc[:, [level_] + samples_list]
    profile_.to_csv(profile_tsv, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser("metagenomics species abundance profiler")
    parser.add_argument("-l", "--abundance_list", type=str, help="abundance list")
    parser.add_argument(
        "--method", default="hsx", choices=["hsx", "jgi"], help="compute method"
    )
    parser.add_argument(
        "--database", default=None, help="contig and genome relationships"
    )
    parser.add_argument(
        "--taxonomy", default=None, help="genome database taxonomy information"
    )
    parser.add_argument("--threads", default=8, type=int, help="threads")
    parser.add_argument("--out_prefix", default="./", type=str, help="output prefix")

    args = parser.parse_args()

    abun_files = pd.read_csv(args.abundance_list, names=["path"]).loc[:, "path"].values

    if args.method == "jgi" and args.database is None:
        print("pleas supply database when parse jgi depth file")
        sys.exit(1)

    if args.method == "hsx":
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, "hsx")
    elif args.method == "jgi":
        profiler_init(args.database)
        count_df, abun_df = get_all_abun_df(abun_files, args.threads, "jgi")
    else:
        print("unsupport method: %s" % args.method)

    outdir = os.path.dirname(args.out_prefix)
    if (outdir == "") or (outdir == "."):
        outdir = "."
    else:
        os.makedirs(outdir, exist_ok=True)
    outprefix = os.path.basename(args.out_prefix)
    count_profile = os.path.join(outdir, outprefix + "_count_profile.tsv")
    abun_profile = os.path.join(outdir, outprefix + "_abundance_profile.tsv")
    count_df.to_csv(count_profile, sep="\t", index=False)
    abun_df.to_csv(abun_profile, sep="\t", index=False)

    if args.taxonomy is not None:
        taxonomy_df = pd.read_csv(args.taxonomy, sep="\t")
        abun_tax_df = abun_df.merge(taxonomy_df)
        samples_list = sorted(abun_df.columns[1:].to_list())

        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_superkingdom_new",
            os.path.join(outdir, outprefix + "_abundance_profile_superkingdom.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_phylum_new",
            os.path.join(outdir, outprefix + "_abundance_profile_phylum.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_order_new",
            os.path.join(outdir, outprefix + "_abundance_profile_order.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_class_new",
            os.path.join(outdir, outprefix + "_abundance_profile_class.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_family_new",
            os.path.join(outdir, outprefix + "_abundance_profile_family.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_genus_new",
            os.path.join(outdir, outprefix + "_abundance_profile_genus.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_species_new",
            os.path.join(outdir, outprefix + "_abundance_profile_species.tsv"),
        )
        get_profile(
            abun_tax_df,
            samples_list,
            "lineages_strain_new",
            os.path.join(outdir, outprefix + "_abundance_profile_strain.tsv"),
        )


if __name__ == "__main__":
    main()
