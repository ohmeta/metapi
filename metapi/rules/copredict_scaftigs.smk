PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                 "gbk", "gff", "sqn", "tbl", "tsv", "txt"]


rule copredict_scaftigs_gene_prodigal:
    input:
        os.path.join(config["output"]["coassembly"],
                     "scaftigs/all.{assembler_co}.out/all.{assembler_co}.scaftigs.fa.gz")
    output:
        pep = os.path.join(
            config["output"]["copredict"],
            "scaftigs_gene/all.{assembler_co}.prodigal.out/all.{assembler_co}.faa"),
        cds = os.path.join(
            config["output"]["copredict"],
            "scaftigs_gene/all.{assembler_co}.prodigal.out/all.{assembler_co}.ffn"),
        gff = os.path.join(
            config["output"]["copredict"],
            "scaftigs_gene/all.{assembler_co}.prodigal.out/all.{assembler_co}.gff")
    log:
        os.path.join(config["output"]["copredict"],
                     "logs/scaftigs_gene/all.{assembler_co}.prodigal.log")
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


rule copredict_scaftigs_gene_prodigal_all:
    input:
        expand(os.path.join(
            config["output"]["copredict"],
            "scaftigs_gene/all.{assembler_co}.prodigal.out/all.{assembler_co}.{ext}"),
               ext=["faa", "ffn", "gff"],
               assembler_co=ASSEMBLERS_CO),

        rules.coassembly_all.input


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
    rule copredict_scaftigs_gene_prokka:
        input:
            os.path.join(
                config["output"]["coassembly"],
                "scaftigs/all.{assembler_co}.out/all.{assembler_co}.scaftigs.fa.gz")
        output:
            expand(os.path.join(
                config["output"]["copredict"],
                "scaftigs_gene/all.{{assembler_co}}.prokka.out/all.{{assembler_co}}.{ext}"),
                   ext=PROKKA_SUFFIX)
        log:
            os.path.join(config["output"]["copredict"],
                         "logs/scaftigs_gene/all.{assembler_co}.prokka.log")
        params:
            prefix = "all.{assembler_co}",
            output_dir = os.path.join(
                config["output"]["copredict"],
                "scaftigs_gene/all.{assembler_co}.prokka.out")
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


    rule copredict_scaftigs_gene_prokka_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["copredict"],
                    "scaftigs_gene/all.{{assembler_co}}.prokka.out/all.{{assembler_co}}.{ext}"),
                ext=PROKKA_SUFFIX)
        output:
            html = os.path.join(
                config["output"]["copredict"],
                "report/scaftigs_gene_{assembler_co}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["copredict"],
                "report/scaftigs_gene_{assembler_co}.multiqc.out/prokka_multiqc_report_data"))
        log:
            os.path.join(
                config["output"]["copredict"],
                "logs/report/scaftigs_gene_{assembler_co}.multiqc.prokka.log")
        params:
            output_dir = os.path.join(
                config["output"]["copredict"],
                "report/scaftigs_gene_{assembler_co}.multiqc.out")
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


    rule copredict_scaftigs_gene_prokka_all:
        input:
            expand([
                os.path.join(
                    config["output"]["copredict"],
                    "scaftigs_gene/all.{assembler_co}.prokka.out/all.{assembler_co}.{ext}"),
                os.path.join(
                    config["output"]["copredict"],
                    "report/scaftigs_gene_{assembler_co}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["copredict"],
                    "report/scaftigs_gene_{assembler_co}.multiqc.out/prokka_multiqc_report_data")],
                   ext=PROKKA_SUFFIX,
                   assembler_co=ASSEMBLERS_CO),

            rules.coassembly_all.input

else:
    rule copredict_scaftigs_gene_prokka_all:
        input:


rule copredict_scaftigs_gene_all:
    input:
        rules.copredict_scaftigs_gene_prodigal_all.input,
        rules.copredict_scaftigs_gene_prokka_all.input


rule predict_scaftigs_gene_all:
    input:
        rules.single_predict_scaftigs_gene_all.input,
        rules.copredict_scaftigs_gene_all.input
