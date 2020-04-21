#!/usr/bin/env python

import argparse
import os
from metapi import configer


def update_config(
    workdir,
    rmhost_host_fasta,
    rmhost_bwa_index,
    rmhost_bowtie2_index,
    kraken2_db,
    prof_index_metadata,
    prof_taxonomy,
    prof_jgi_index,
):

    conf_file = os.path.join(workdir, "config.yaml")
    conf_file_up = os.path.join(workdir, "config_update.yaml")

    conf = configer.parse_yaml(os.path.join(workdir, "config.yaml"))

    conf["params"]["rmhost"]["host_fasta"] = rmhost_host_fasta
    conf["params"]["rmhost"]["bwa"]["index_prefix"] = rmhost_bwa_index
    conf["params"]["rmhost"]["bowtie2"]["index_prefix"] = rmhost_bowtie2_index
    conf["params"]["classify"]["kraken2"]["database"] = kraken2_db
    conf["params"]["profiling"]["index_metadata"] = prof_index_metadata
    conf["params"]["profiling"]["taxonomy"] = prof_taxonomy
    conf["params"]["profiling"]["jgi"]["index_prefix"] = prof_jgi_index

    configer.update_config(conf_file, conf_file_up, conf, remove=False)
    os.rename(conf_file_up, conf_file)


def main():
    parser = argparse.ArgumentParser("update metapi config.yaml")
    parser.add_argument("-d", "--workdir", type=str, help="work dir", default="./")
    parser.add_argument("-a", "--rmhost_host_fasta", type=str, help="rmhost host fasta")
    parser.add_argument(
        "-i", "--rmhost_bwa_index", type=str, help="rmhost bwa index prefix"
    )
    parser.add_argument(
        "-I", "--rmhost_bowtie2_index", type=str, help="rmhost bowtie2 index prefix"
    )
    parser.add_argument("-k", "--kraken2_db", type=str, help="kraken2 database")
    parser.add_argument(
        "-m", "--profiling_index_metadata", type=str, help="profiling index metadata"
    )
    parser.add_argument(
        "-t", "--profiling_taxonomy", type=str, help="profiling taxonomy"
    )
    parser.add_argument(
        "-j", "--profiling_jgi_index", type=str, help="profiling jgi index prefix"
    )
    args = parser.parse_args()

    update_config(
        args.workdir,
        args.rmhost_host_fasta,
        args.rmhost_bwa_index,
        args.rmhost_bowtie2_index,
        args.kraken2_db,
        args.profiling_index_metadata,
        args.profiling_taxonomy,
        args.profiling_jgi_index,
    )


if __name__ == "__main__":
    main()
