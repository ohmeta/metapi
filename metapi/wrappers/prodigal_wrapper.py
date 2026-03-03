#!/usr/bin/env python3

import glob
import os
import stat
import sys
import subprocess
import concurrent.futures

import pandas as pd
from checkm import prodigal


def _validate_prodigal_output(pep_file, cds_file, gff_file, bin_id):
    """Validate Prodigal output files exist and are non-empty."""
    for f, fname in [(pep_file, "pep"), (cds_file, "cds"), (gff_file, "gff")]:
        if not os.path.exists(f):
            print(f"Warning: Prodigal did not generate {fname} file for {bin_id}")
            return False
        if os.path.getsize(f) == 0:
            print(f"Warning: Prodigal generated empty {fname} file for {bin_id}")
            return False
    return True


def _cleanup_output_files(pep_file, cds_file, gff_file, best_translation_table=None):
    """Remove output files from failed Prodigal runs."""
    files_to_clean = []
    if best_translation_table in [4, 11]:
        files_to_clean = [
            f"{pep_file}.{best_translation_table}",
            f"{cds_file}.{best_translation_table}",
            f"{gff_file}.{best_translation_table}"
        ]
    else:
        files_to_clean = [pep_file, cds_file, gff_file]
    for f in files_to_clean:
        if os.path.exists(f):
            subprocess.run(f'''rm -rf {f}''', shell=True)


def _infer_translation_table(pep_file, cds_file, gff_file):
    """Infer translation table from existing output files with .4 or .11 suffix."""
    for tt in [4, 11]:
        if (os.path.exists(f"{pep_file}.{tt}") or 
            os.path.exists(f"{cds_file}.{tt}") or 
            os.path.exists(f"{gff_file}.{tt}")):
            return tt
    return None


def _compress_files(pep_file, cds_file, gff_file):
    """Compress output files using pigz."""
    for f in [pep_file, cds_file, gff_file]:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            subprocess.run(f'''pigz -f {f}''', shell=True)


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

    if not os.path.exists(bin_fa):
        print(f"Warning: Input file does not exist for {bin_id}: {bin_fa}")
        return None

    prodigal_runner = prodigal.ProdigalRunner(output_dir)
    prodigal_runner.aaGeneFile = pep_file
    prodigal_runner.ntGeneFile = cds_file
    prodigal_runner.gffFile = gff_file

    best_translation_table = None
    try:
        best_translation_table = prodigal_runner.run(bin_fa, True)
    except SystemExit as e:
        print(f"Warning: Prodigal exited with error for {bin_id}: {e}")
        inferred_tt = _infer_translation_table(pep_file, cds_file, gff_file)
        _cleanup_output_files(pep_file, cds_file, gff_file, inferred_tt)
        return None
    except Exception as e:
        print(f"Warning: Prodigal failed for {bin_id} with exception: {e}")
        inferred_tt = _infer_translation_table(pep_file, cds_file, gff_file)
        _cleanup_output_files(pep_file, cds_file, gff_file, inferred_tt)
        return None

    if best_translation_table is None:
        print(f"Warning: Prodigal returned no translation table for {bin_id}")
        _cleanup_output_files(pep_file, cds_file, gff_file, None)
        return None

    if not _validate_prodigal_output(pep_file, cds_file, gff_file, bin_id):
        print(f"Warning: Prodigal failed or returned no output for {bin_id} (code: 0 but no valid output)")
        _cleanup_output_files(pep_file, cds_file, gff_file, best_translation_table)
        return None

    _compress_files(pep_file, cds_file, gff_file)

    if not os.path.exists(pep_file_gz) or os.stat(pep_file_gz)[stat.ST_SIZE] == 0:
        print(f"Warning: Compressed file is missing or empty for {bin_id}")
        _cleanup_output_files(pep_file_gz, cds_file_gz, gff_file_gz, None)
        return None

    if best_translation_table in [4, 11]:
        return (bin_id, bin_fa, pep_file_gz, best_translation_table)
    else:
        return (bin_id, bin_fa, pep_file_gz, f"unknown: {best_translation_table}")


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
    for result in executor.map(run_prodigal, input_list):
        if result is not None:
            table_list.append(result)

print("Prodigal annotation completed for all bins.")
table_df = pd.DataFrame(table_list, columns=["bin_id", "bin_file", "pep_file", "best_translation_table"])
table_df.to_csv(output_done, sep="\t", index=False)