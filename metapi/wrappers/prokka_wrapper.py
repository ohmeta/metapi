import glob
import os
import time
import subprocess


PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                 "gbk", "gff", "sqn", "tbl", "tsv", "txt"]

bin_list = glob.glob(snakemake.input["mags_dir"] + "/*.fa")
gff_count = 0

for bin_fa in bin_list:
    bin_id = os.path.basename(os.path.splitext(bin_fa)[0])
    output_dir = os.path.join(snakemake.params["output_dir"], bin_id)
    gff_file = os.path.join(output_dir, bin_id + ".gff")

    subprocess(f'''echo "\nProcessing {bin_fa}\n" >> {snakemake.log}''', shell=True)
    subprocess(
        f'''
        prokka {bin_fa} \
        --force \
        --centre X \
        --compliant \
        --cpus {snakemake.threads} \
        --outdir {output_dir} \
        --locustag {bin_id} \
        --prefix {bin_id} \
        --kingdom {snakemake.params["kingdom"]} \
        2>> {snakemake.log} 
        ''', shell=True)

    if os.path.exists(gff_file):
        gff_count += 1

if gff_count == len(bin_list):
    subprocess('''touch {snakemake.output["done"]}''', shell=True)

    for suffix in PROKKA_SUFFIX:
        prokka_f = os.path.join(output_dir, f'''{bin_id}.{suffix}''')
        if os.path.exists(prokka_f):
            subprocess.run(f'''pigz {prokka_f}''', shell=True)