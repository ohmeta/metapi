if config["params"]["classify"]["kraken2"]["do"]:
    rule classify_short_reads_kraken2:
        input:
            assembly_input
        output:
            table = protected(os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.table.gz")),
            report = protected(os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.report.gz"))
        log:
            os.path.join(config["output"]["classify"],
                         "logs/{sample}.kraken2.log")
        params:
            paired = "--paired" if IS_PE else "",
            database = config["params"]["classify"]["kraken2"]["database"],
            use_mpa_style = "--use-mpa-style" \
                if config["params"]["classify"]["kraken2"]["use_mpa_style"] \
                   else "",
            report_zero_counts = "--report-zero-counts" \
                if config["params"]["classify"]["kraken2"]["report_zero_counts"] \
                   else ""
        threads:
            config["params"]["classify"]["threads"]
        run:
            shell(
                '''
                kraken2 \
                --use-names \
                --threads {threads} \
                --db {params.database} \
                --output {output.table} \
                --report {output.report} \
                {params.use_mpa_style} \
                {params.report_zero_counts} \
                {params.paired} \
                --gzip-compressed \
                {input} \
                2> {log}
                ''')

            table = os.path.splitext(output.table)[0]
            report = os.path.splitext(output.table)[0]

            if os.path.exists(table) and os.path.exists(report):
                shell('''pigz -p {threads} %s''' \
                      % os.path.splitext(output.table)[0])
                shell('''pigz -p {threads} %s''' \
                      % os.path.splitext(output.report)[0])


    rule classify_short_reads_kraken2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.table.gz"),
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.report.gz")],
                   sample=SAMPLES.index.unique()),

            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule classify_short_reads_kraken2_all:
        input:


if config["params"]["classify"]["gtdbtk"]["do"]:
    rule classify_hmq_bins_gtdbtk:
        input:
            os.path.join(
                config["output"]["checkm"],
                "bins_hmq/{assembler}.{binner}.links")
        output:
            directory(os.path.join(
                config["output"]["classify"],
                "bins_hmq/{assembler}.{binner}.gtdbtk.out"))
        log:
            os.path.join(config["output"]["classify"],
                         "logs/{assembler}.{binner}.gtdbtk.log")
        params:
            bin_suffix = "fa"
        threads:
            config["params"]["classify"]["threads"]
        shell:
            '''
            gtdbtk classify_wf \
            --genome_dir {input}/ \
            --out_dir {output} \
            --extension {params.bin_suffix} \
            --cpus {threads} \
            > {log}
            '''


    rule classify_hmq_bins_gtdbtk_all:
        input:
            expand(
                os.path.join(
                    config["output"]["classify"],
                    "bins_hmq/{assembler}.{binner}.gtdbtk.out"),
                assembler=ASSEMBLERS,
                binner=BINNERS),

            rules.checkm_all.input,

else:
    rule classify_hmq_bins_gtdbtk_all:
        input:


rule classify_all:
    input:
        rules.classify_short_reads_kraken2_all.input,
        rules.classify_hmq_bins_gtdbtk_all.input
