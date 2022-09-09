import os
import pandas as pd
import subprocess
import concurrent.futures


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
gtdb_tables = snakemake.input["gtdb_table"]
rep_genomes_info = snakemake.input["rep_genomes_info"]
gtdb_to_ncbi_script = snakemake.params["gtdb_to_ncbi_script"]
metadata_archaea = snakemake.params["metadata_archaea"]
metadata_bacteria = snakemake.params["metadata_bacteria"]
table_gtdb = snakemake.output["table_gtdb"]
table_ncbi = snakemake.output["table_ncbi"]
table_all = snakemake.output["table_all"]


gtdb_list = []
ncbi_list = []

for i in gtdb_tables:
    out_dir = os.path.dirname(i)
    ar122_tsv = os.path.join(out_dir, "gtdbtk.ar122.summary.tsv")
    bac120_tsv = os.path.join(out_dir, "gtdbtk.bac120.summary.tsv")

    if os.path.exists(ar122_tsv):
        gtdb_list.append(ar122_tsv)
    if os.path.exists(bac120_tsv):
        gtdb_list.append(bac120_tsv)

    gtdb_to_ncbi_summary = os.path.join(out_dir, "gtdbtk.ncbi.summary.tsv")
    gtdb_to_ncbi_log = os.path.join(out_dir, "gtdbtk.to.ncbi.log")

    gtdb_to_ncbi_cmd = \
        f'''
        python {gtdb_to_ncbi_script} \
        --gtdbtk_output_dir {out_dir} \
        --output_file {gtdb_to_ncbi_summary} \
        --ar122_metadata_file {metadata_archaea} \
        --bac120_metadata_file {metadata_bacteria} \
        > {gtdb_to_ncbi_log}
        '''

    if os.path.exists(gtdb_to_ncbi_summary):
        ncbi_list.append(gtdb_to_ncbi_summary)

merge(gtdb_list, parse, threads, output=table_gtdb)
merge(ncbi_list, parse, threads, output=table_ncbi)

table_gtdb = pd.read_csv(table_gtdb, sep="\t").rename(columns={"classification": "GTDB classification"})
table_ncbi = pd.read_csv(table_ncbi, sep="\t")
table_rep_genomes_info = pd.read_csv(rep_genomes_info, sep="\t").rename(columns={"genome": "user_genome"})

table_all = pd.merge(table_gtdb, table_ncbi, how="inner", on=["user_genome", "GTDB classification"]).\
                merge(table_rep_genomes_info, how="inner", on="user_genome")

table_all.to_csv(table_all, sep="\t", index=False)
