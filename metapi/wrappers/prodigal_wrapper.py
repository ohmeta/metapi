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

    bin_id = os.path.basename(os.path.splitext(os.path.splitext(bin_fa)[0])[0])

    pep_file = os.path.join(output_dir, bin_id + ".faa")
    cds_file = os.path.join(output_dir, bin_id + ".ffn")
    gff_file = os.path.join(output_dir, bin_id + ".gff")

    pep_file_gz = pep_file + ".gz"
    cds_file_gz = cds_file + ".gz"
    gff_file_gz = gff_file + ".gz"

    prodigal_runner = prodigal.ProdigalRunner(output_dir)
    prodigal_runner.aaGeneFile = pep_file
    prodigal_runner.ntGeneFile = cds_file
    prodigal_runner.gffFile = gff_file

    best_translation_table = prodigal_runner.run(bin_fa, True)

    if os.path.exists(pep_file) and (os.path.getsize(pep_file) > 0):
        subprocess.run(f'''pigz -f {pep_file}''', shell=True)
        if os.path.exists(cds_file) and (os.path.getsize(cds_file) > 0):
            subprocess.run(f'''pigz -f {cds_file}''', shell=True)
        if os.path.exists(gff_file) and (os.path.getsize(gff_file) > 0):
            subprocess.run(f'''pigz -f {gff_file}''', shell=True)
    else:
        subprocess.run(f'''rm -rf {pep_file}''', shell=True)
        subprocess.run(f'''rm -rf {cds_file}''', shell=True)
        subprocess.run(f'''rm -rf {gff_file}''', shell=True)
 
    if best_translation_table in [4, 11]:
        if (os.path.exists(pep_file_gz)) and (os.path.exists(cds_file_gz)) and (os.path.exists(gff_file_gz)) and (os.stat(pep_file_gz)[stat.ST_SIZE]) > 0:
            return (bin_id, bin_fa, pep_file_gz, best_translation_table)
        else:
            return None
    else:
        if (os.path.exists(pep_file_gz)) and (os.path.exists(cds_file_gz)) and (os.path.exists(gff_file_gz)) and (os.stat(pep_file_gz)[stat.ST_SIZE]) > 0:
            return (bin_id, bin_fa, pep_file_gz, f"unknown: {best_translation_table}")
        else: 
            return None


workers = int(sys.argv[1])
input_mags_dir = os.path.dirname(sys.argv[2])
output_done = sys.argv[3]
output_dir = os.path.dirname(output_done)

bin_list = glob.glob(input_mags_dir + "/*.fa.gz")

input_list = []
for bin_fa in bin_list:
    input_list.append((bin_fa, output_dir))

table_list = []


subprocess.run(f'''rm -rf {output_dir}''', shell=True)
subprocess.run(f'''mkdir -p {output_dir}''', shell=True)

with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
    for table_df in executor.map(run_prodigal, input_list):
        if table_df is not None:
            table_list.append(table_df)

table_df = pd.DataFrame(table_list, columns=["bin_id", "bin_file", "pep_file", "best_translation_table"])
table_df.to_csv(output_done, sep="\t", index=False)