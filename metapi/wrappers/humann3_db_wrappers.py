#!/usr/bin/env python3

import argparse
import logging
import time
import os

from humann import config
from humann.humann import timestamp_message
from humann.search import prescreen
from humann.search import nucleotide


def build_chocophlan_pangenome_db(
    log_file, file_basename, db_dir, prescreen_threshold, taxonomic_profile
):
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=log_file,
        format="%(asctime)s - %(name)s - %(levelname)s: %(message)s",
        level=getattr(logging, config.log_level),
        filemode="w",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    config.file_basename = file_basename
    config.temp_dir = db_dir
    config.unnamed_temp_dir = db_dir
    config.prescreen_threshold = prescreen_threshold

    start_time = time.time()
    custom_database = prescreen.create_custom_database(
        config.nucleotide_database, taxonomic_profile
    )
    start_time = timestamp_message(
        "parse taxonomy profile and custom database creation", start_time
    )

    if custom_database != "Empty":
        nucleotide_index_file = nucleotide.index(custom_database)
    start_time = timestamp_message("database index", start_time)

    os.remove(custom_database)


def main():
    parser = argparse.ArgumentParser(
        description="HUMANn3 build chocophlan pangenome database wrapper"
    )
    parser.add_argument("--log", help="log file")
    parser.add_argument("--basename", type=str, help="file basename")
    parser.add_argument("--db_dir", type=str, help="database directory")
    parser.add_argument("--prescreen_threshold", type=float, help="prescreen threshold")
    parser.add_argument(
        "--taxonomic_profile", type=str, help="MetaPhlAn3 taxonomic profile"
    )

    args = parser.parse_args()
    build_chocophlan_pangenome_db(
        args.log,
        args.basename,
        args.db_dir,
        args.prescreen_threshold,
        args.taxonomic_profile,
    )


if __name__ == "__main__":
    main()
