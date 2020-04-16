#!/usr/bin/env python3

import os
import pandas as pd
import ncbi_genome_download as ngd


def download_genomes(config):
    """
    ncbi-genome-download
    """
    ngd.download(
        format="fasta,assembly-report",
        assembly_level="complete",
        taxid=config["params"]["simulation"]["taxid"],
        refseq_category="reference",
        output_folder=config["results"]["simulation"]["genomes"],
        retries=3,
        metadata_table=os.path.join(
            config["results"]["simulation"]["genomes"], "metadata.tsv"
        ),
        group="bacteria",
    )


def generate_samples(config):
    """
    generate samples
    """
    download_genomes(config)

    samples_list = []
    for i in config["params"]["simulation"]["samples"]:
        sample_id = i[0]
        fq_prefix = "_".join(i)
        fq1 = os.path.join(
            config["results"]["simulation"]["short_reads"], fq_prefix + ".1.fq.gz"
        )
        fq2 = os.path.join(
            config["results"]["simulation"]["short_reads"], fq_prefix + ".2.fq.gz"
        )
        samples_list.append([sample_id, fq1, fq2])

    samples_df = pd.DataFrame(samples_list, columns=["id", "fq1", "fq2"]).set_index(
        "id", drop=False
    )
    return samples_df
