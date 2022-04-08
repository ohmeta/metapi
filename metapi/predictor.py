#!/usr/bin/env python

import glob
import os
import stat
import time
import subprocess
from Bio import SeqIO


def calculate_coding_density(gff):
    pass


def run_prodigal(bin_fa, mode, condon_table, log_file, pep_file, cds_file, gff_file):
    pep = f'''{pep_file}.{condon_table}'''
    cds = f'''{cds_file}.{condon_table}'''
    gff = f'''{gff_file}.{condon_table}'''

    cmd = f'''prodigal -i {bin_fa} -m -a {pep} -d {cds} -o {gff} -f gff -g {condon_table} -p {mode} 2>> {log_file}'''

    subprocess.run('''echo "Predict g ene from {bin_fa} using condon translation {condon_table}\n >> {log_file}''')
    subprocess.run(cmd, shell=True)

    if (not os.path.exists(pep)) or (os.stat(pep)[stat.ST_SIZE] == 0):
        if mode == "single":
            cmd = cmd.replace("-p single", "-p meta")
            subprocess.run(cmd, shell=True)

    if (os.path.exists(pep)) and (os.stat(pep)[stat.ST_SIZE] > 0):
        return True
    else:
        return False


def predict_bins(bins_list, output_dir, log_file):

    for bin_fa in bins_list:
        bin_id = os.path.basename(os.path.splitext(bin_fa)[0])

        pep_file = os.path.join(output_dir, bin_id + ".faa")
        cds_file = os.path.join(output_dir, bin_id + ".ffn")
        gff_file = os.path.join(output_dir, bin_id + ".gff")
        
        coding_density = {4: 0, 11: 0}
        total_bases = 0
        mode = "single"
        for seq in SeqIO.parse(bin_fa, "fasta"):
            total_bases += len(seq)
        if total_bases < 100000:
            mode = "meta"

        for condon_table in [4, 11]:
            done = run_prodigal(bin_fa, mode, condon_table, log_file, pep_file, cds_file, gff_file)
            if done:
                gff = f'''{gff_file}.{condon_table}'''
                coding_density_info = calculate_coding_density(gff)