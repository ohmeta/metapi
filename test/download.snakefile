dwn_lnks = {
    '1': 'https://molb7621.github.io/workshop/_downloads/sample.fa',
    '2': 'https://molb7621.github.io/workshop/_downloads/sample.fa'
}

import os

# association between chromosomes and their links
def chromo2link(wildcards):
    return dwn_lnks[wildcards.chromo]

rule all:
    input:
        os.path.join('data/genome_dir', 'human_en37_sm.fa')

rule download:
    output:
        "data/chr_dir/{chromo}"
    params:
        link=chromo2link
    shell:
        "wget {params.link} -O {output}"

rule merger:
    input:
        expand(os.path.join('data/chr_dir', "{chromo}"),
               chromo=dwn_lnks.keys())
    output:
        os.path.join('data/genome_dir', 'human_en37_sm.fa')
    run:
        txt = open("{output}", 'a+')
        with open(os.path.join('data/chr_dir', "{chromo}"), 'r') as file:
            line = file.readline()
            while line:
                txt.write(line)
                line = file.read(line)
        txt.close()
