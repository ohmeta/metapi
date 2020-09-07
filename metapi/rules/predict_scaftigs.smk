PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                 "gbk", "gff", "sqn", "tbl", "tsv", "txt"]


rule predict_scaftigs_gene_prodigal:
    input:
        os.path.join(config["output"]["assembly"],
                     "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        pep = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.faa"),
        cds = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.ffn"),
        gff = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.gff")
    log:
        os.path.join(config["output"]["predict"],
                     "logs/scaftigs_gene/{sample}.{assembler}.prodigal.log")
    params:
        format = config["params"]["predict"]["format"]
    shell:
        '''
        zcat {input} | \
        prodigal \
        -m \
        -a {output.pep} \
        -d {output.cds} \
        -o {output.gff} \
        -f {params.format} \
        -p meta -q \
        2> {log}
        '''


rule predict_scaftigs_gene_prodigal_all:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.{ext}"),
               ext=["faa", "ffn", "gff"],
               assembler=ASSEMBLERS,
               sample=SAMPLES.index.unique()),

        rules.single_assembly_all.input


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
    rule predict_scaftigs_gene_prokka:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{sample}}.{{assembler}}.prokka.out/{{sample}}.{{assembler}}.{ext}"),
                   ext=PROKKA_SUFFIX)
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{sample}.{assembler}.prokka.log")
        params:
            prefix = "{sample}.{assembler}",
            output_dir = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prokka.out")
        threads:
            config["params"]["predict"]["threads"]
        shell:
            '''
            prokka {input} \
            --force \
            --centre X \
            --compliant \
            --cpus {threads} \
            --outdir {params.output_dir} \
            --locustag {params.prefix} \
            --prefix {params.prefix} \
            --metagenome \
            2> {log}
            '''


    rule predict_scaftigs_gene_prokka_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{{assembler}}.prokka.out/{sample}.{{assembler}}.{ext}"),
                ext=PROKKA_SUFFIX,
                sample=SAMPLES.index.unique())
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data"))
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/report/scaftigs_gene_{assembler}.multiqc.prokka.log")
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out")
        shell:
            '''
            multiqc \
            --cl_config "prokka_fn_snames: True" \
            --outdir {params.output_dir} \
            --title prokka \
            --module prokka \
            {input} \
            2> {log}
            '''


    rule predict_scaftigs_gene_prokka_all:
        input:
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{assembler}.prokka.out/{sample}.{assembler}.{ext}"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data")],
                   ext=PROKKA_SUFFIX,
                   assembler=ASSEMBLERS,
                   sample=SAMPLES.index.unique()),

            rules.single_assembly_all.input

else:
    rule predict_scaftigs_gene_prokka_all:
        input:


rule single_predict_scaftigs_gene_all:
    input:
        rules.predict_scaftigs_gene_prodigal_all.input,
        rules.predict_scaftigs_gene_prokka_all.input
