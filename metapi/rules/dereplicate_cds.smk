rule dereplicate_gene_prepare:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{sample}.{{assembler}}.prodigal.out/{sample}.{{assembler}}.ffn"),
               sample=SAMPLES.index.unique())
    output:
        os.path.join(config["output"]["predict"],
                     "{assembler}.prodigal.scaftigs.gene.merged.ffn")
    shell:
        '''
        cat {input} > {output}
        '''


if config["params"]["dereplicate"]["cdhit"]["do_gene"]:
    rule dereplicate_gene_cdhit:
        input:
            os.path.join(config["output"]["predict"],
                         "{assembler}.prodigal.scaftigs.gene.merged.ffn")
        output:
            os.path.join(config["output"]["dereplicate"],
                         "genes/{assembler}.prodigal.scaftigs.gene.merged.nr.ffn")
        log:
            os.path.join(config["output"]["dereplicate"],
                         "logs/{assembler}.prodigal.scaftigs.gene.cdhit.log")
        threads:
            config["params"]["dereplicate"]["cdhit"]["threads"]
        params:
            identity = config["params"]["dereplicate"]["cdhit"]["identity"],
            overlap = config["params"]["dereplicate"]["cdhit"]["overlap"],
            wordlen = config["params"]["dereplicate"]["cdhit"]["wordlen"],
            global_ = 1 if config["params"]["dereplicate"]["cdhit"]["global"] else 0,
            memory = config["params"]["dereplicate"]["cdhit"]["memory"],
            clstrlen = config["params"]["dereplicate"]["cdhit"]["clstrlen"],
            default_algorithm = 0 if config["params"]["dereplicate"]["cdhit"]["default_algorithm"] else 1,
            both_alignment = 1 if config["params"]["dereplicate"]["cdhit"]["both_alignment"] \
                else 0
        shell:
            '''
            cd-hit-est -i {input} -o {output} \
            -c {params.identity} \
            -n {params.wordlen} \
            -G {params.global_} \
            -aS {params.overlap} \
            -M {params.memory} \
            -d {params.clstrlen} \
            -g {params.default_algorithm} \
            -r {params.both_alignment} \
            -T {threads} 2> {log}
            '''


    rule dereplicate_gene_cdhit_all:
        input:
            expand(os.path.join(
                config["output"]["dereplicate"],
                "genes/{assembler}.prodigal.scaftigs.gene.merged.nr.ffn"),
                   assembler=ASSEMBLERS)

else:
    rule dereplicate_gene_cdhit_all:
        input:


rule dereplicate_gene_all:
    input:
        rules.dereplicate_gene_cdhit_all.input
