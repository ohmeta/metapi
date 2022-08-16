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


if IS_PE:
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
                f'''
                metaspades.py \
                --continue \
                -o {snakemake.params["output_dir"]} \
                > {snakemake.log}
                ''', shell=True)
    else:
        subprocess(f'''rm -rf {snakemake.params["output_dir"]}''', shell=True)
        subprocess(f'''mkdir -p {snakemake.params["output_dir"]}.reads''', shell=True)

        reads_num = len(snakemake.input)
        if reads_num == 2:
            subprocess(
                f'''
                metaspades.py \
                -1 {snakemake.input[0]} \
                -2 {snakemake.input[1]} \
                -k {snakemake.params["kmers"]} \
                {snakemake.params["only_assembler"]} \
                --memory {snakemake.params["memory"]} \
                --threads {snakemake.threads} \
                --checkpoints last \
                -o {snakemake.params["output_dir"]} \
                > {snakemake.log}
                ''', shell=True)
        else:
            subprocess(f'''cat {snakemake.input[0:reads_num//2]} > {snakemake.params["output_dir"]}.reads/reads.1.fq.gz''',
                       shell=True)
            subprocess(f'''cat {snakemake.input[reads_num//2:]} > {snakemake.params["output_dir"]}.reads/reads.2.fq.gz''',
                       shell=True)

            subprocess(
                f'''
                metaspades.py \
                -1 {snakemake.params["output_dir"]}.reads/reads.1.fq.gz \
                -2 {snakemake.params["output_dir"]}.reads/reads.2.fq.gz \
                -k {snakemake.params["kmers"]} \
                {snakemake.params["only_assembler"]} \
                --memory {snakemake.params["memory"]} \
                --threads {snakemake.threads} \
                --checkpoints last \
                -o {snakemake.params["output_dir"]} \
                > {snakemake.log}
                ''', shell=True)

            subprocess('''rm -rf {snakemake.params["output_dir"]}.reads''', shell=True)

    subprocess(
        '''
        pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/scaffolds.fasta && \
        mv {snakemake.params["output_dir"]}/scaffolds.fasta.gz \
        {snakemake.params["output_dir"]}/{snakemake.params["prefix"].metaspades.scaffolds.fa.gz
        ''', shell=True)
    subprocess(
        '''
        pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/contigs.fasta && \
        mv {snakemake.params["output_dir"]}/contigs.fasta.gz \
        {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.metaspades.contigs.fa.gz
        ''', shell=True)
    subprocess(
        '''
        pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/contigs.paths && \
        mv {snakemake.params["output_dir"]}/contigs.paths.gz \
        {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.metaspades.contigs.paths.gz
        ''', shell=True)
    subprocess(
        '''
        pigz -p {snakemake.threads} {snakemake.params["output_dir"]}/scaffolds.paths && \
        mv {snakemake.params["output_dir"]}/scaffolds.paths.gz \
        {snakemake.params["output_dir"]}/{snakemake.params["prefix"]}.metaspades.scaffolds.paths.gz
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
            ln -s {snakemake.params["prefix"]}.metaspades.scaffolds.fa.gz \
            {snakemake.params["prefix"]}.metaspades.scaftigs.fa.gz && \
            ln -s {snakemake.params["prefix"]}.metaspades.scaffolds.paths.gz \
            {snakemake.params["prefix"]}.metaspades.scaftigs.paths.gz && \
            popd
            ''', shell=True)
    else:
            subprocess(
            '''
            pushd {snakemake.params["output_dir"]} && \
            ln -s {snakemake.params["prefix"]}.metaspades.contigs.fa.gz \
            {snakemake.params["prefix"]}.metaspades.scaftigs.fa.gz && \
            ln -s {snakemake.params["prefix"]}.metaspades.contigs.paths.gz \
            {snakemake.params["prefix"]}.metaspades.scaftigs.paths.gz && \
            popd
            ''', shell=True)

    if snakemake.params["only_save_scaftigs"]:
        subprocess('''fd -d 1 -E "*.gz" . {snakemake.params["output_dir"]} -x rm -rf {{}}''', shell=True)
    else:
        subprocess(
            '''
            rm -rf {snakemake.params["output_dir"]}/{corrected, misc, pipeline_state, tmp} && \
            tar -czvf {snakemake.params["tar_results"]}.gz {snakemake.params["output_dir"]}/K*
            ''', shell=True)
else:
    print(
        '''
        Don't support single-end reads assembly using MetaSPAdes\n,
        you can try SPAdes or MegaHit, IDBA_UD
        ''')
    sys.exit(1)
