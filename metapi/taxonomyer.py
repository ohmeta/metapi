#!/usr/bin/env python

import os
import gzip
import pandas as pd
from Bio import SeqIO, bgzf
import pandas as pd


def gen_genome_index(count, max_num=6):
    count_str = str(count)
    count_str_list = [x for x in count_str]
    if len(count_str_list) <= max_num:
        count_str_list = ["0"] * (max_num - len(count_str_list)) + count_str_list
    return "".join(count_str_list)


def set_genomes(rep_info, map_name, base_dir):
    genome_range = len(str(len(rep_info)))
    if genome_range <= 3:
        max_num = 3
    elif genome_range <= 6:
        max_num = 6
    elif genome_range <= 9:
        max_num = 9
    elif genome_range <= 12:
        max_num = 12
    else:
        max_num = 15

    for i in range(0, len(rep_info)):
        genome_index = gen_genome_index(i, max_num=max_num)

        out_dir = os.path.realpath(base_dir)
        for k in range(0, max_num, 3):
            out_dir = os.path.join(out_dir, genome_index[k:k+3])

        genome_id = f"{map_name}_GENOME_{genome_index}"
        genome_path = os.path.join(out_dir, f"{genome_id}.fna.gz")

        rep_info.at[i, "genome_id"] = genome_id
        rep_info.at[i, "genome_path"] = genome_path
    return rep_info


def set_lineages(genome_id, classification, rep_level):
    LINEAGES = ["k", "p", "c", "o", "f", "g", "s", "t"]

    tax_list = []
    for i in classification.split(";"):
        tax_list.append(i.split("__", maxsplit=1)[1].replace(" ", "_"))
    ##print(tax_list)

    lineage_list = []

    if len(tax_list) == 7:
        lineage_list = tax_list[0:6]
        if rep_level == "strain":
            lineage_list.append(tax_list[6])
        elif rep_level == "species":
            lineage_list.append(f"{tax_list[6]}-{genome_id}")
    elif len(tax_list) < 7:
        lineage_list = tax_list
        for tax_index in range(len(tax_list), 7):
            tax_level = LINEAGES[tax_index-1]
            if tax_index == 6:
                if rep_level == "strain":
                    lineage_list.append(f"{tax_level}__{lineage_list[-1]}-unclassified")
                elif rep_level == "species":
                    lineage_list.append(f"{tax_level}__{lineage_list[-1]}-{genome_id}")
            else:
                lineage_list.append(f"{tax_level}__{lineage_list[-1]}-unclassified")

    lineage_list.append(f"s__{lineage_list[-1]}-{genome_id}")
    ##print(lineage_list)

    lineage = ";".join(lineage_list)
    ##print(lineage)

    lineage_list.append(lineage)
    ##print(f"####### length of lineages list: {len(lineage_list)}")

    return pd.Series(lineage_list)


def update_genomes(rep_info):
    for i in range(0, len(rep_info)):
        clade_lineage = rep_info.at[i, "lineage"]
        completeness = rep_info.at[i, "completeness"]
        contamination = rep_info.at[i, "contamination"]
        strain_heterogeneity = rep_info.at[i, "strain_heterogeneity"]
        quality_score = rep_info.at[i, "quality_score"]
        bin_file = rep_info.at[i, "bin_file"]
        bin_id = rep_info.at[i, "bin_id"]
        genome_path = rep_info.at[i, "genome_path"]
        genome_id = rep_info.at[i, "genome_id"]

        if bin_file.endswith(".gz"):
            handle = gzip.open(bin_file, 'rt')
        else:
            handle = open(bin_file, 'r')

        os.makedirs(os.path.dirname(genome_path), exist_ok=True)

        with bgzf.BgzfWriter(genome_path, 'wb') as oh:
            j = -1
            for rc in SeqIO.parse(handle, "fasta"):
                j += 1
                contig_name = f"{genome_id}_{j}"

                rc_id = rc.id
                rc.id = contig_name
                rc.description = f"{genome_id}|original_contig_id={rc_id}|original_bin_id={bin_id}|gtdb_classification={clade_lineage}|completeness={completeness}|contamination={contamination}|strain_heterogeneity={strain_heterogeneity}|quality_scroe={quality_score}"
                SeqIO.write(rc, oh, "fasta") 
 
        handle.close()


def refine_taxonomy(genomes_info_f, tax_info_f, map_name, rep_level, base_dir, out_file):
    genomes_info = pd.read_csv(genomes_info_f, sep="\t")
    tax_info = pd.read_csv(tax_info_f, sep="\t").rename(columns={"user_genome": "genome"})

    rep_info = pd.merge(genomes_info, tax_info, how="inner", on=["genome"])\
                 .sort_values(["classification", "genome"])\
                 .reset_index(drop=True)

    # step 1: set genomes
    rep_info = set_genomes(rep_info, map_name, base_dir)

    # step 2: set lineages
    rep_info[["kingdom", "phylum", "class", "order",
              "family", "genus", "species", "strain",
              "lineage"]] = rep_info.apply(
                lambda x: set_lineages(
                    x["genome_id"], x["classification"], rep_level),
                    axis = 1, result_type = "expand")

    # step 3: save taxonomy
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    rep_info.to_csv(out_file, sep="\t", index=False)

    # setp 4: update genomes
    update_genomes(rep_info)

