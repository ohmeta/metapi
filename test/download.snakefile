dwn_lnks = {
    '1': 'https://molb7621.github.io/workshop/_downloads/sample.fa',
    '2': 'https://molb7621.github.io/workshop/_downloads/sample.fa'
}

import os

def get_link(wildcards):
    return dwn_lnks[wildcards.chromo]

rule all:
    input:
        os.path.join('data/genome_dir', 'human_en37_sm.fa')

rule download:
    output:
        os.path.join('data/chr_dir', "{chromo}")
    params:
        link = lambda wildcards: dwn_lnks[wildcards.chromo]
    shell:
        "wget {params.link} -O {output}"

rule merger:
    input:
        expand(os.path.join('data/chr_dir', "{chromo}"), chromo=dwn_lnks.keys())
    output:
        os.path.join('data/genome_dir', 'human_en37_sm.fa')
    run:
        txt = open(output[0], "a+")
        for i in input:
            with open(i, "r") as file:
                line = file.readline()
                while line:
                    txt.write(line)
                    line = file.readline()
            txt.write("\n")
        txt.close()
