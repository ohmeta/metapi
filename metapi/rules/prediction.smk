
if config["params"]["prediction"]["prodigal"]["do"]:
    rule prediction:
        input:
            scaftigs = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            pep = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.pep.faa"),
            cds = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.cds.ffn"),
            gff = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.cds.gff"),
            start = os.path.join(config["results"]["prediction"], "{sample}.prodigal_out/{sample}.score.gff")
        log:
            os.path.join(config["logs"]["prediction"], "{sample}.prodigal.log")
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
            -s {output.start} \
            -f {params.format} -p {params.mode} -q \
            2> {log}
            '''