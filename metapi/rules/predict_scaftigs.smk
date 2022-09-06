rule predict_scaftigs_gene_prodigal:
    input:
        os.path.join(config["output"]["assembly"],
                     "scaftigs/{binning_group}.{assembly_group}.{assembler}.out/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        multiext(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal.out/{binning_group}.{assembly_group}.{assembler}"),
            "faa", "ffn", "gff")
    conda:
        config["envs"]["predict"]
    log:
        os.path.join(config["output"]["predict"],
                     "logs/scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal.log")
    shell:
        '''
        zcat {input} | \
        prodigal \
        -m \
        -a {output[0]} \
        -d {output[1]} \
        -o {output[2]} \
        -f gff \
        -p meta -q \
        2> {log}
        '''


rule predict_scaftigs_gene_prodigal_all:
    input:
        expand(expand(
            os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal.out/{binning_group}.{assembly_group}.{assembler}.{{ext}}"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
                ext = ["faa", "ffn", "gff"]),

        rules.assembly_all.input


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
    PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                     "gbk", "gff", "sqn", "tbl", "tsv", "txt"]


    rule predict_scaftigs_gene_prokka:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}.out/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{binning_group}}.{{assembly_group}}.{{assembler}}.prokka.out/{{binning_group}}.{{assembly_group}}.{{assembler}}.{ext}"),
                ext=PROKKA_SUFFIX)
        conda:
            config["envs"]["predict"]
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka.log")
        params:
            prefix = "{binning_group}.{assembly_group}.{assembler}",
            output_dir = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka.out")
        threads:
            config["params"]["predict"]["threads"]
        shell:
            '''
            gzip -dc {input} > {params.output_dir}/input.fa

            prokka {params.output_dir}/input.fa \
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
            expand(expand(
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{binning_group}.{assembly_group}.{{{{assembler}}}}.prokka.out/{binning_group}.{assembly_group}.{{{{assembler}}}}.{{ext}}"),
                    zip,
                    binning_group=ASSEMBLY_GROUP["binning_group"],
                    assembly_group=ASSEMBLY_GROUP["assembly_group"]),
                    ext=PROKKA_SUFFIX)
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data"))
        conda:
            config["envs"]["multiqc"]
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
            expand(expand(
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka.out/{binning_group}.{assembly_group}.{assembler}.{{ext}}"),
                    zip,
                    binning_group=ASSEMBLY_GROUPS["binning_group"],
                    assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                    assembler=ASSEMBLY_GROUPS["assembler"]),
                    ext = PROKKA_SUFFIX),
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data")],
                   assembler=ASSEMBLERS),

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