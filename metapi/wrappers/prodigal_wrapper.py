#!/usr/bin/env python3

import glob
import os
import stat
import sys
import subprocess
import concurrent.futures

import pandas as pd
from checkm import prodigal


def run_prodigal(input_list):
    bin_fa = os.path.abspath(input_list[0])
    output_dir = os.path.abspath(input_list[1])

    bin_id = os.path.basename(os.path.splitext(bin_fa)[0])

    pep_file = os.path.join(output_dir, bin_id + ".faa")
    cds_file = os.path.join(output_dir, bin_id + ".ffn")
    gff_file = os.path.join(output_dir, bin_id + ".gff")

    prodigal_runner = prodigal.ProdigalRunner(output_dir)
    prodigal_runner.aaGeneFile = pep_file
    prodigal_runner.ntGeneFile = cds_file
    prodigal_runner.gffFile = gff_file

    best_translation_table = prodigal_runner.run(bin_fa, True)

    if (best_translation_table in [4, 11]) and (os.path.exists(pep_file)) and (os.stat(pep_file)[stat.ST_SIZE]) > 0:
        return (bin_id, bin_fa, pep_file, best_translation_table)
    else:
        return (bin_id, bin_fa, pep_file, f"unknown: {best_translation_table}")


workers = int(sys.argv[1])
input_bins_dir = os.path.dirname(sys.argv[2])
output_done = sys.argv[3]
output_dir = os.path.dirname(output_done)

bin_list = glob.glob(input_bins_dir + "/*.fa")

input_list = []
for bin_fa in bin_list:
    input_list.append((bin_fa, output_dir))

table_list = []


subprocess.run(f'''rm -rf {output_dir}''', shell=True)
subprocess.run(f'''mkdir -p {output_dir}''', shell=True)

with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
    for table_df in executor.map(run_prodigal, input_list):
        table_list.append(table_df)

table_df = pd.DataFrame(table_list, columns=["bin_id", "bin_file", "pep_file", "best_translation_table"])
table_df.to_csv(output_done, sep="\t", index=False)