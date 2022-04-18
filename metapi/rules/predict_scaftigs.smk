PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                 "gbk", "gff", "sqn", "tbl", "tsv", "txt"]


rule predict_scaftigs_gene_prodigal:
    input:
        os.path.join(config["output"]["assembly"],
                     "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        pep = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.faa"),
        cds = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.ffn"),
        gff = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.gff")
    log:
        os.path.join(config["output"]["predict"],
                     "logs/scaftigs_gene/{assembly_group}.{assembler}.prodigal.log")
    shell:
        '''
        zcat {input} | \
        prodigal \
        -m \
        -a {output.pep} \
        -d {output.cds} \
        -o {output.gff} \
        -f gff \
        -p meta -q \
        2> {log}
        '''


rule predict_scaftigs_gene_prodigal_all:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.{ext}"),
               ext=["faa", "ffn", "gff"],
               assembler=ASSEMBLERS,
               assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST),

        rules.assembly_all.input


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
    rule predict_scaftigs_gene_prokka:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{assembly_group}}.{{assembler}}.prokka.out/{{assembly_group}}.{{assembler}}.{ext}"),
                   ext=PROKKA_SUFFIX)
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{assembly_group}.{assembler}.prokka.log")
        params:
            prefix = "{assembly_group}.{assembler}",
            output_dir = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{assembly_group}.{assembler}.prokka.out")
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
                    "scaftigs_gene/{assembly_group}.{{assembler}}.prokka.out/{assembly_group}.{{assembler}}.{ext}"),
                ext=PROKKA_SUFFIX,
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
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
                    "scaftigs_gene/{assembly_group}.{assembler}.prokka.out/{assembly_group}.{assembler}.{ext}"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data")],
                   ext=PROKKA_SUFFIX,
                   assembler=ASSEMBLERS,
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST),

            rules.assembly_all.input

else:
    rule predict_scaftigs_gene_prokka_all:
        input:


rule predict_scaftigs_gene_all:
    input:
        rules.predict_scaftigs_gene_prodigal_all.input,
        rules.predict_scaftigs_gene_prokka_all.input


localrules:
    predict_scaftigs_gene_prodigal_all,
    predict_scaftigs_gene_prokka_all,
    predict_scaftigs_gene_all,