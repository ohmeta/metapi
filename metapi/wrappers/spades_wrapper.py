
import os
import re
import sys
import subprocess


def parse_spades_params(params_file):
    with open(params_file, "r") as ih:
        cmd = ih.readline().strip()

        matches = re.match(r".*-k\t(.*?)\t--memory\t(\d+)\t--threads\t(\d+).*", cmd)
        if matches:
            kmers = str(matches.group(1))
            memory = str(matches.group(2))
            threads = str(matches.group(3))
            if "--only-assembler" in cmd:
                return [kmers, memory, threads, True]
            else:
                return [kmers, memory, threads, False]
        else:
            return None


continue_assembly = False
error_correct = False if snakemake.params["only_assembler"] == "" else True
params_file = os.path.join(snakemake.params["output_dir"], "params.txt")

if os.path.exists(params_file):
    spades_params = parse_spades_params(params_file)
    if spades_params is not None:
        if all([snakemake.params["kmers"] == spades_params[0],
                snakemake.params["memory"] == spades_params[1],
                str(snakemake.threads) == spades_params[2],
                error_correct == spades_params[3]]):
            continue_assembly = True


if continue_assembly:
    subprocess(
        '''
        spades.py \
        --continue \
        -o {snakemake.params["output_dir"]} \
        > {snakemake.log}
        ''', shell=True)
else:
    subprocess('''rm -rf {snakemake.params["output_dir"]}''', shell=True)
    reads_num = len(snakemake.input)
    input_str = ""
    if IS_PE:
        if reads_num == 2:
            input_str = f'''-1 {snakemake.input[0]} -2 {snakemake.input[1]}'''
        else:
            subprocess('''makedir -p {snakemake.params["output_dir"]}.reads''')
            subprocess(f'''cat {snakemake.input[0:reads_num//2]} > {snakemake.params["output_dir"]}.reads/reads.1.fq.gz''')
            subprocess(f'''cat {snakemake.input[reads_num//2:]} > {snakemake.params["output_dir"]}.reads/reads.2.fq.gz''')
            input_str = f'''-1 {snakemake.params["output_dir"]}.reads/reads.1.fq.gz -2 {snakemake.params["output_dir"]}.reads/reads.2.fq.gz'''
    else:
        if reads_num == 1:
            input_str = f'''-s {snakemake.input[0]}'''
        else:
            subprocess('''makedir -p {snakemake.params["output_dir"]}.reads''')
            subprocess('''cat {input} > {snakemake.params["output_dir"]}.reads/reads.fq.gz''')
            input_str = f'''-s {snakemake.params["output_dir"]}.reads/reads.fq.gz'''

    subprocess(
        f'''
        spades.py \
        {input_str} \
        -k {snakemake.params["kmers"]} \
        {snakemake.params["only_assembler"]} \
        --memory {snakemake.params["memory"]} \
        --threads {snakemake.threads} \
        --checkpoints last \
        -o {snakemake.params["output_dir"]} \
        > {snakemake.log}
        ''', shell=True)

subprocess(
    '''
    pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/scaffolds.fasta && \
    mv {snakemake.params["output_dir"]}/scaffolds.fasta.gz \
    {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.spades.scaffolds.fa.gz
    ''', shell=True)
subprocess(
    '''
    pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/contigs.fasta && \
    mv {snakemake.params["output_dir"]}/contigs.fasta.gz \
    {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.spades.contigs.fa.gz
    ''', shell=True)
subprocess(
    '''
    pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/contigs.paths && \
    mv {snakemake.params["output_dir"]}/contigs.paths.gz \
    {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.spades.contigs.paths.gz
    ''', shell=True)
subprocess(
    '''
    pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/scaffolds.paths && \
    mv {snakemake.params["output_dir"]}/scaffolds.paths.gz \
    {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.spades.scaffolds.paths.gz
    ''', shell=True)
subprocess(
    '''
    pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/assembly_graph_with_scaffolds.gfa && \
    mv {snakemake.params["output_dir"]}/assembly_graph_with_scaffolds.gfa.gz {snakemake.output["gfa"]}
    ''', shell=True)

if snakemake.params["link_scaffolds"]:
    subprocess(
        '''
        pushd {snakemake.params["output_dir"]} && \
        ln -s {snakemake.params["prefix"]}.spades.scaffolds.fa.gz \
        {snakemake.params["prefix"]}.spades.scaftigs.fa.gz && \
        ln -s {snakemake.params["prefix"]}.spades.scaffolds.paths.gz \
        {snakemake.params["prefix"]}.spades.scaftigs.paths.gz && \
        popd
        ''', shell=True)
else:
    subprocess(
        '''
        pushd {snakemake.params["output_dir"]} && \
        ln -s {snakemake.params["prefix"]}.spades.contigs.fa.gz \
        {snakemake.params["prefix"]}.spades.scaftigs.fa.gz && \
        ln -s {snakemake.params["prefix"]}.spades.contigs.paths.gz \
        {snakemake.params["prefix"]}.spades.scaftigs.paths.gz && \
        popd
        ''', shell=True)

if snakemake.params["only_save_scaftigs"]:
    subprocess('''fd -d 1 -E "*.gz" . {snakemake.params["output_dir"]} -x rm -rf {{}}''', shell=True)
else:
    subprocess(
        '''
        rm -rf {snakemake.params["output_dir"]}/{corrected, misc, pipeline_state, tmp} &&
        tar -czvf {snakemake.params["tar_results"]}.gz {snakemake.params["output_dir"]}/K*
        ''', shell=True)
