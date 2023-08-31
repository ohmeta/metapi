import os
import pandas as pd
import subprocess
import concurrent.futures
from pprint import pprint


def parse(stats_file):
    if os.path.exists(stats_file):
        try:
            df = pd.read_csv(stats_file, sep="\t")
        except pd.errors.EmptyDataError:
            print("%s is empty, please check" % stats_file)
            return None

        if not df.empty:
            return df
        else:
            return None
    else:
        print("%s is not exists" % stats_file)
        return None


def merge(input_list, func, workers, **kwargs):
    df_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        for df in executor.map(func, input_list):
            if df is not None:
                df_list.append(df)

    df_ = pd.concat(df_list)

    if "output" in kwargs:
        df_.to_csv(kwargs["output"], sep="\t", index=False)
    return df_


threads = int(snakemake.threads)

gtdb_done_list = snakemake.input["gtdb_done"]

gtdb_to_ncbi_script = snakemake.params["gtdb_to_ncbi_script"]
metadata_archaea = snakemake.params["metadata_archaea"]
metadata_bacteria = snakemake.params["metadata_bacteria"]

table_gtdb = snakemake.output["table_gtdb"]
table_ncbi = snakemake.output["table_ncbi"]
table_all = snakemake.output["table_all"]

os.makedirs(os.path.dirname(table_all), exist_ok=True)

gtdb_list = []
ncbi_list = []

for i in gtdb_done_list:
    out_dir = os.path.dirname(i)
    archaea_tsv = os.path.join(out_dir, "gtdbtk.archaea.summary.tsv")
    bacteria_tsv = os.path.join(out_dir, "gtdbtk.bacteria.summary.tsv")

    if os.path.exists(archaea_tsv):
        gtdb_list.append(archaea_tsv)
    if os.path.exists(bacteria_tsv):
        gtdb_list.append(bacteria_tsv)

    gtdb_to_ncbi_summary = os.path.join(out_dir, "gtdbtk.ncbi.summary.tsv")
    gtdb_to_ncbi_log = os.path.join(out_dir, "gtdbtk.to.ncbi.log")

    archaea_cmd = "--ar53_metadata_file"
    if "ar122" in os.path.realpath(archaea_tsv):
        archaea_cmd = "--ar122_metadata_file"

    bacteria_cmd = "--bac120_metadata_file"

    gtdb_to_ncbi_cmd = \
        f'''
        python {gtdb_to_ncbi_script} \
        --gtdbtk_output_dir {out_dir} \
        --output_file {gtdb_to_ncbi_summary} \
        {archaea_cmd} {metadata_archaea} \
        {bacteria_cmd} {metadata_bacteria} \
        > {gtdb_to_ncbi_log}
        '''
    subprocess.run(gtdb_to_ncbi_cmd, shell=True)

    if os.path.exists(gtdb_to_ncbi_summary):
        ncbi_list.append(gtdb_to_ncbi_summary)


if len(gtdb_list) > 0:
    table_gtdb_df = merge(gtdb_list, parse, threads, output=table_gtdb)
else:
    print(f"No {table_gtdb} generate")

if len(ncbi_list) > 0:
    table_ncbi_df = merge(ncbi_list, parse, threads, output=table_ncbi)
else:
    print(f"No {table_ncbi} generate")


table_gtdb_df = table_gtdb_df.rename(columns={"classification": "GTDB classification"})
pprint(table_gtdb_df)

table_ncbi_df = table_ncbi_df.rename(columns={"Genome ID": "user_genome"})
pprint(table_ncbi_df)

table_all_df = pd.merge(
    table_gtdb_df, table_ncbi_df, how="inner",
    on=["user_genome", "GTDB classification"])#\

table_all_df.to_csv(table_all, sep="\t", index=False)
