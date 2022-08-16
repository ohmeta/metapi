#!/usr/bin/env python

import subprocess
import os


subprocess('''rm -rf {snakemake.params["output_dir"]}''', shell=True)
subprocess('''mkdir {snakemake.params["output_dir"]}''', shell=True)

reads = os.path.join(
    snakemake.params["output_dir"],
    f'''{snakemake.params["prefix"]}.fa''')

reads_num = len(snakemake.input)

if IS_PE:
    if reads_num == 2:
        subprocess(
            f'''
            seqtk mergepe {snakemake.input[0]} {snakemake.input[1]} | \
            seqtk seq -A - > {reads}
            ''', shell=True)
    else:
        subprocess(
            f'''
            cat {snakemake.input[0:reads_num//2]} > {reads}.1.fq.gz
            cat {snakemake.input[reads_num//2:]} > {reads}.2.fq.gz
            seqtk mergepe {reads}.1.fq.gz {reads}.2.fq.gz | \
            seqtk seq -A - > {reads}
            rm -rf {reads}.1.fq.gz {reads}.2.fq.gz
            ''', shell=True)
else:
    if reads_num == 1:
        subprocess(f'''seqtk seq -A {snakemake.input[0]} > {reads}''', shell=True)
    else:
        subprocess(
            f'''
            cat {snakemake.input} > {reads}.fq.gz
            seqtk seq -A {snakemake.input[0]} > {reads}
            rm -rf {reads}.fq.gz
            ''', shell=True)

subprocess(
    f'''
    idba_ud \
    -r {reads} \
    --mink {snakemake.params["mink"]} \
    --maxk {snakemake.params["maxk"]} \
    --step {snakemake.params["step"]} \
    --min_contig {snakemake.params["min_contig"]} \
    -o {snakemake.params["output_dir"]} \
    --num_snakemake.threads {snakemake.threads} \
    --pre_correction \
    > {snakemake.log}
    ''', shell=True)

subprocess(f'''rm -rf {reads}''', shell=True)
subprocess(f'''sed -i 's#^>#>{snakemake.params["prefix"]}_#g' {snakemake.params["output_dir"]}/scaffold.fa''', shell=True)
subprocess(f'''pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/scaffold.fa''', shell=True)
subprocess(f'''mv {snakemake.params["output_dir"]}/scaffold.fa.gz {snakemake.output["scaftigs"]}''', shell=True)

if snakemake.params["only_save_scaftigs"]:
    subprocess(
        '''
        find {snakemake.params["output_dir"]} \
        -type f \
        ! -wholename "{snakemake.output['scaftigs']}" -delete
        ''', shell=True)
