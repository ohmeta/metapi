rule predict_scaftigs_gene_prodigal:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal.{ext}.gz"),
            ext=["faa", "ffn", "gff"])
    log:
        os.path.join(
            config["output"]["predict"],
            "logs/predict_scaftigs_gene_prodigal/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["predict"],
            "benchmark/predict_scaftigs_gene_prodigal/{binning_group}.{assembly_group}.{assembler}.txt")
    conda:
        config["envs"]["predict"]
    threads:
        config["params"]["predict"]["threads"]
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

        pigz -f -p {threads} ${{FAA%.gz}}
        pigz -f -p {threads} ${{FFN%.gz}}
        pigz -f -p {threads} ${{GFF%.gz}}
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


PROKKA_SUFFIX = [
    "err", "log", "faa", "ffn", "fna", "fsa",
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
    log:
        os.path.join(
            config["output"]["predict"],
            "logs/predict_scaftigs_gene_prokka/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["predict"],
            "benchmark/predict_scaftigs_gene_prokka/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        prefix = "{binning_group}.{assembly_group}.{assembler}.prokka",
        output_dir = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prokka")
    threads:
        config["params"]["predict"]["threads"]
    conda:
        config["envs"]["predict"]
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

        pigz -f -p {threads} ${{ERR%.gz}}
        pigz -f -p {threads} ${{LOG%.gz}}
        pigz -f -p {threads} ${{FAA%.gz}}
        pigz -f -p {threads} ${{FFN%.gz}}
        pigz -f -p {threads} ${{FNA%.gz}}
        pigz -f -p {threads} ${{FSA%.gz}}
        pigz -f -p {threads} ${{GBK%.gz}}
        pigz -f -p {threads} ${{GFF%.gz}}
        pigz -f -p {threads} ${{SQN%.gz}}
        pigz -f -p {threads} ${{TBL%.gz}}
        pigz -f -p {threads} ${{TSV%.gz}}
        pigz -f -p {threads} ${{TXT%.gz}}
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
    log:
        os.path.join(
            config["output"]["predict"],
            "logs/predict_scaftigs_gene_prokka_multiqc/scaftigs_gene_{assembler}.multiqc.prokka.log")
    benchmark:
        os.path.join(
            config["output"]["predict"],
            "benchmark/predict_scaftigs_gene_prokka_multiqc/scaftigs_gene_{assembler}.multiqc.prokka.log")
    params:
        output_dir = os.path.join(
            config["output"]["predict"],
            "report/scaftigs_gene_{assembler}.multiqc")
    conda:
        config["envs"]["multiqc"]
    shell:
        '''
        multiqc \
        --cl_config "prokka_fn_snames: True" \
        --outdir {params.output_dir} \
        --title prokka \
        --module prokka \
        {input} \
        >{log} 2>&1
        '''


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
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