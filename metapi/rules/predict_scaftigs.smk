rule predict_scaftigs_gene_prodigal:
    input:
        os.path.join(config["output"]["assembly"],
                     "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal.{ext}.gz"),
            ext=["faa", "ffn", "gff"])
    conda:
        config["envs"]["predict"]
    benchmark:
        os.path.join(config["output"]["predict"],
                     "benchmark/scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal.benchmark.txt")
    log:
        os.path.join(config["output"]["predict"],
                     "logs/scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal.log")
    shell:
        '''
        FAA={output[0]}
        FFN={output[1]}
        GFF={output[2]}

        zcat {input} | \
        prodigal \
        -m \
        -a ${{FAA%.gz}} \
        -d ${{FFN%.gz}} \
        -o ${{GFF%.gz}} \
        -f gff \
        -p meta \
        -q \
        >{log} 2>&1

        pigz -f ${{FAA%.gz}}
        pigz -f ${{FFN%.gz}}
        pigz -f ${{GFF%.gz}}
        '''


rule predict_scaftigs_gene_prodigal_all:
    input:
        expand(expand(
            os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binning_group}.{assembly_group}.{assembler}.prodigal.{{ext}}.gz"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
                ext = ["faa", "ffn", "gff"])


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
    PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                     "gbk", "gff", "sqn", "tbl", "tsv", "txt"]


    rule predict_scaftigs_gene_prokka:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{binning_group}}.{{assembly_group}}.{{assembler}}.prokka/{{binning_group}}.{{assembly_group}}.{{assembler}}.prokka.{ext}.gz"),
                ext=PROKKA_SUFFIX)
        conda:
            config["envs"]["predict"]
        benchmark:
            os.path.join(config["output"]["predict"],
                         "benchmark/scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka.benchmark.txt")
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka.log")
        params:
            prefix = "{binning_group}.{assembly_group}.{assembler}.prokka",
            output_dir = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka")
        threads:
            config["params"]["predict"]["threads"]
        shell:
            '''
            ERR={output[0]}
            LOG={output[1]}
            FAA={output[2]}
            FFN={output[3]}
            FNA={output[4]}
            FSA={output[5]}
            GBK={output[6]}
            GFF={output[7]}
            SQN={output[8]}
            TBL={output[9]}
            TSV={output[10]}
            TXT={output[11]}

            prokka \
            <(zcat {input}) \
            --force \
            --centre X \
            --compliant \
            --cpus {threads} \
            --outdir {params.output_dir} \
            --locustag {params.prefix} \
            --prefix {params.prefix} \
            --metagenome \
            2> {log}

            pigz -f ${{ERR%.gz}}
            pigz -f ${{LOG%.gz}}
            pigz -f ${{FAA%.gz}}
            pigz -f ${{FFN%.gz}}
            pigz -f ${{FNA%.gz}}
            pigz -f ${{FSA%.gz}}
            pigz -f ${{GBK%.gz}}
            pigz -f ${{GFF%.gz}}
            pigz -f ${{SQN%.gz}}
            pigz -f ${{TBL%.gz}}
            pigz -f ${{TSV%.gz}}
            pigz -f ${{TXT%.gz}}
            '''


    rule predict_scaftigs_gene_prokka_multiqc:
        input:
            expand(expand(
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{binning_group}.{assembly_group}.{{{{assembler}}}}.prokka/{binning_group}.{assembly_group}.{{{{assembler}}}}.prokka.{{ext}}.gz"),
                    zip,
                    binning_group=ASSEMBLY_GROUP["binning_group"],
                    assembly_group=ASSEMBLY_GROUP["assembly_group"]),
                    ext=PROKKA_SUFFIX)
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc/prokka_multiqc_report_data"))
        conda:
            config["envs"]["multiqc"]
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/report/scaftigs_gene_{assembler}.multiqc.prokka.log")
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc")
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
                    "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{binning_group}.{assembly_group}.{assembler}.prokka.{{ext}}.gz"),
                    zip,
                    binning_group=ASSEMBLY_GROUPS["binning_group"],
                    assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                    assembler=ASSEMBLY_GROUPS["assembler"]),
                    ext = PROKKA_SUFFIX),
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc/prokka_multiqc_report_data")],
                   assembler=ASSEMBLERS)

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
    predict_scaftigs_gene_all