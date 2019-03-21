
if config["params"]["prediction"]["prodigal"]["do"]:
    rule prediction:
        input:
            os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            pep = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.{assembler}.pep.faa"),
            cds = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.{assembler}.cds.ffn"),
            gff = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.{assembler}.cds.gff"),
            score = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.{assembler}.score.gff")
        log:
            os.path.join(config["logs"]["prediction"], "{sample}.{assembler}.prodigal.log")
        params:
            format = config["params"]["prediction"]["prodigal"]["format"],
            mode = config["params"]["prediction"]["prodigal"]["mode"]
        shell:
            '''
            prodigal \
            -i {input} \
            -a {output.pep} \
            -d {output.cds} \
            -o {output.gff} \
            -s {output.score} \
            -f {params.format} -p {params.mode} -q \
            2> {log}
            '''