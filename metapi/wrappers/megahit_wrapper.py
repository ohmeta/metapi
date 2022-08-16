#!/usr/bin/env python

from Bio import SeqIO
import re
import os
import subprocess


input_str = ""
kmer_opts = ""
reads_num = 0 

if snakemake.params["presets"] != "":
    kmer_opts = f'''--presets {snakemake.params["presets"]}'''
else:
    kmer_opts = f'''--k-list {snakemake.params["k_list"]}'''

if os.path.exists(os.path.join(snakemake.params["output_dir"], "options.json")):
    subprocess.run(f'''megahit --continue --out-dir {snakemake.params["output_dir"]}''',
                   shell=True)
else:
    reads_num = len(snakemake.input)

if IS_PE:
    if reads_num == 2:
        input_str = f'''-1 {snakemake.input[0]} -2 {snakemake.input[1]}'''
    else:
        input_str = f'''-1 {",".join(snakemake.input[0:reads_num//2])} -2 {",".join(snakemake.input[reads_num//2:])}'''
else:
    input_str = f'''-r {",".join(snakemake.input)}'''

subprocess.run(f'''rm -rf {snakemake.params["output_dir"]}''', shell=True)

subprocess.run(
    f'''
    megahit \
    {input_str} \
    -t {snakemake.threads} \
    {kmer_opts} \
    --min-contig-len {snakemake.params["min_contig"]} \
    --out-dir {snakemake.params["output_dir"]} \
    --out-prefix {snakemake.params["output_prefix"]} \
    2> {snakemake.log}
    ''', shell=True)

k_num = 0
for seq_record in SeqIO.parse(snakemake.params["contigs"], "fasta"):
    k_num = int(re.search('k(.*)_', seq_record.id).group(1))
    break

subprocess.run(
    f'''
    megahit_toolkit contig2fastg \
    {k_num} \
    {snakemake.params["contigs"]} \
    > {snakemake.params["fastg"]}
    ''', shell=True)

subprocess.run(f'''fastg2gfa {snakmeake.params["fastg"]} > {snakemake.params["gfa"]}''', shell=True)
subprocess.run(f'''pigz -p {snakemake.threads} {snakemake.params["fastg"]}''', shell=True)
subprocess.run(f'''pigz -p {snakemake.threads} {snakemake.params["gfa"]}''', shell=True)

subprocess.run(f'''pigz -p {snakemake.threads} {snakemake.params["contigs"]}''', shell=True)
subprocess.run(f'''mv {snakemake.params["contigs"]}.gz {snakemake.output["scaftigs"]}''', shell=True)

if snakemake.params["only_save_scaftigs"]:
    subprocess.run(f'''fd -t f -E "*.gz" . {snakemake.params["output_dir"]} -x rm -rf {{}}''', shell=Tru)
    subprocess.run(f'''rm -rf {snakemake.params["output_dir"]}/intermediate_contigs''', sehll=True)
