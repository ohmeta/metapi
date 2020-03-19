#!/usr/bin/env python

import argparse
import bz2
import gzip
import os
import pickle
import sys

import magic
import ncbi_genome_download as ngd
import pandas as pd


def set_kingdom(x):
    return x["lineages"].split("|")[0].split("__")[-1].lower()


def set_genus(x):
    return x["lineages"].split("|")[5].split("__")[-1].replace("_", " ")


def set_species(x):
    return x["lineages"].split("|")[-1].split("__")[-1].replace("_", " ")


def parse_pickle(mpa_pickle):
    magic_str = magic.from_file(mpa_pickle)
    if "gzip" in magic_str:
        handle = gzip.open(mpa_pickle, "rb")
    elif "bzip2" in magic_str:
        handle = bz2.open(mpa_pickle, "rb")
    elif "data" in magic_str:
        handle = open(mpa_pickle, "rb")
    else:
        print(
            """%s: unknown file type: gzip? bzip2? pickle?\n"""
            """Please check it first""" % mpa_pickle
        )
        sys.exit(1)

    mpa_dict = {}
    mpa = pickle.load(handle)

    for marker in mpa["markers"]:
        tax_id = marker.split("__")[0]
        if tax_id not in mpa_dict:
            mpa_dict[tax_id] = mpa["markers"][marker]["taxon"]
    mpa_df = pd.DataFrame(mpa_dict.items(), columns=["taxid", "lineages"]).set_index(
        "lineages"
    )
    return mpa_df


def parse_taxonomy(taxonomy):
    tax = pd.read_csv(taxonomy, sep="\t").iloc[:, 0].to_frame(name="lineages")
    tax["kingdom"] = tax.apply(lambda x: set_kingdom(x), axis=1)
    tax["genus"] = tax.apply(lambda x: set_genus(x), axis=1)
    tax["species"] = tax.apply(lambda x: set_species(x), axis=1)
    return tax.set_index("lineages")


def main():
    parser = argparse.ArgumentParser(description="ncbi-genome-downloader wrapper")
    parser.add_argument(
        "--mpa_pkl", help="MetaPhlAn2 clade pickle",
    )
    parser.add_argument("--tax_list", help="a file contain species taxonmy lineages")
    parser.add_argument("--print", action="store_true", help="print download shell")
    parser.add_argument("--download", action="store_true", help="just download")
    parser.add_argument("--outdir", default=None, type=str, help="output directory")
    parser.add_argument(
        "--parallel", default=4, type=int, help="parallel download number, default: 4"
    )
    parser.add_argument(
        "--logdir", default="./ngd_logs", type=str, help="ngd download log directory"
    )
    args = parser.parse_args()

    mpa_df = parse_pickle(args.mpa_pkl)
    tax_df = parse_taxonomy(args.tax_list)

    if args.print:
        for i in tax_df.index:
            cmd = """ngd \
            --section genbank \
            --format all \
            --assembly-level all
            --genus "{genus}" \
            --species-taxid {taxid} \
            --refseq-category representative \
            --output-folder {outdir} \
            --flat-output \
            --parallel {parallel} \
            --retries 3 \
            --metadata-table {table}
            {group}""".format(
                genus=tax_df.loc[i, "genus"],
                taxid=mpa_df.loc[i, "taxid"],
                outdir=args.outdir,
                parallel=args.parallel,
                table=os.path.join(args.logdir, i + ".ngd.log"),
                group=tax_df.loc[i, "kingdom"],
            )
            print(cmd)

    if args.download:
        os.mkdirs(args.outdir, exists_ok=True)
        os.mkdirs(args.logdir, exists_ok=True)
        for i in tax_df.index:
            ngd.download(
                section="genbank",
                format="all",
                assembly_level="all",
                genus=tax_df.loc[i, "genus"],
                species_taxid=mpa_df.loc[i, "taxid"],
                refseq_category="representative",
                outdir=args.outdir,
                flat_output=True,
                parallel=args.parallel,
                retries=3,
                metadata_table=os.path.join(args.logdir, i + ".ngd.log"),
                group=tax_df.loc[i, "kingdom"],
            )


if __name__ == "__main__":
    main()
