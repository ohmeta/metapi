#!/usr/bin/env python3

import os
import pandas as pd
import concurrent.futures


def gen_samples_info(samples, output, config):
    """
    GSC MIxS environmental sample template
    * fileds are mandatory
    """

    header = [
        "*sample_name",
        "sample_title",
        "*tax_id",
        "*organism",
        "*investigation_type",
        "*collection_date",
        "*env_biome",
        "*env_feature",
        "*env_material",
        "*geo_loc_name",
        "*host",
        "*lat_lon",
        "*strain",
        "*estimated_size",
        "*isol_growth_condt",
        "*num_replicons",
        "*ref_biomaterial",
        "*ploidy",
        "*propagation",
        "biotic_relationship",
        "hem_administration" "encoded_traits",
        "ethnicity",
        "extrachrom_elements",
        "health_state",
        "host_age",
        "host_body_mass_index",
        "host_body_product",
        "host_body_temp",
        "host_diet",
        "host_disease",
        "host_family_relationship",
        "host_genotype",
        "host_height",
        "host_last_meal",
        "host_occupation",
        "host_phenotype",
        "host_pulse",
        "host_sex",
        "host_subject_id",
        "host_taxid",
        "host_tissue_sampled",
        "host_tot_mass",
        "ihmc_medication_code",
        "isolation_source",
        "medic_hist_perform",
        "misc_param",
        "nose_mouth_teeth_throat_disord",
        "organism_count",
        "oxy_stat_samp",
        "pathogenicity",
        "perturbation",
        "rel_to_oxygen",
        "samp_collect_device",
        "samp_mat_process",
        "samp_salinity",
        "samp_size",
        "samp_store_dur",
        "samp_store_loc",
        "samp_store_temp",
        "samp_vol_we_dna_ext",
        "source_material_id",
        "subspecf_gen_lin",
        "temp",
        "time_last_toothbrush",
        "trophic_level",
        "description",
    ]

    samples_df = pd.DataFrame(columns=header)
    samples_df["*sample_name"] = samples.index.unique()
    samples_df["sample_title"] = samples.index.unique()
    samples_df["source_material_id"] = samples.index.unique()
    for key in config["upload"]["samples"]:
        samples_df["*" + key] = config["upload"]["samples"][key]
    samples_df.to_excel(output, index=False, startrow=12)


def parse_md5(md5_file):
    try:
        if os.path.exists(md5_file):
            df = pd.read_csv(
                md5_file, sep="\s+", header=None, names=["file_md5", "file_name"]
            )
            if df.empty:
                print("%s is empty, please check" % md5_file)
                return None
            df["sample_name"] = df.apply(
                lambda x: os.path.basename(x["file_name"]).split(".")[0], axis=1
            )
            df["file_name"] = df.apply(
                lambda x: os.path.basename(x["file_name"]), axis=1
            )
            if len(df) == 2:
                df_fq1 = df.iloc[0].to_frame().T
                df_fq2 = (
                    df.iloc[1]
                    .to_frame()
                    .T.rename(
                        columns={"file_name": "file2_name", "file_md5": "file2_md5"}
                    )
                )
                return pd.merge(df_fq1, df_fq2)
            else:
                df["file2_name"] = ""
                df["file2_md5"] = ""
                return df
        else:
            print("%s is not exists" % md5_file)
            return None
    except pd.io.common.EmptyDataError:
        print("%s is empty, please check" % md5_file)
        return None


def gen_info(input_list, output, config, workers, group):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df_ in executor.map(parse_md5, input_list):
            if df_ is not None:
                df_list.append(df_)
    df = pd.concat(df_list)

    for key in config["upload"][group].keys():
        df[key] = config["upload"][group][key]

    df["project_accession"] = config["upload"]["project_accession"]

    if group == "sequencing_run":
        run_df = df.loc[
            :,
            ["project_accession", "sample_name"]
            + list(config["upload"][group].keys())
            + ["file_name", "file_md5", "file2_name", "file2_md5"],
        ]
        if config["params"]["reads_layout"] == "pe":
            run_df["library_layout"] = "paired"
        else:
            run_df["library_layout"] = "fragment/single"
            run_df["nominal_size"] = ""

        run_df["library_name"] = run_df["sample_name"]

        run_df = run_df.rename(
            columns={
                "project_accession": "*project_accession",
                "sample_name": "*sample_name",
                "experiment_title": "*experiment_title",
                "library_name": "*library_name",
                "library_strategy": "*library_strategy",
                "library_selection": "*library_selection",
                "library_layout": "*library_layout",
                "platform": "*platform",
                "instrument_model": "*instrument_model",
                "spot_layout": "*spot_layout",
                "file_type": "*file_type",
                "file_name": "*file_name",
                "file_md5": "*file_md5",
                "file2_name": "*file2_name",
                "file2_md5": "*file2_md5",
            }
        )
        run_df.to_excel(output, sheet_name="Metadata", index=False)

    if group == "assembly":
        asm_df = df.rename(
            columns={"file_name": "fasta_file_name", "file_md5": "fasta_file_md5"}
        )
        asm_df["assembly_name"] = df["sample_name"] + "-" + df["assembly_method"]
        asm_df["sample_accession"] = ""
        asm_df = asm_df.loc[
            :,
            ["project_accession", "sample_accession", "sample_name", "assembly_name"]
            + list(config["upload"][group].keys())
            + ["fasta_file_name", "fasta_file_md5"],
        ]
        asm_df.to_excel(output, sheet_name="Genome_Assembly", index=False)
