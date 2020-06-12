if config["params"]["classify"]["kraken2"]["do"]:
    rule classify_short_reads_kraken2:
        input:
            reads = assembly_input,
            database = expand(os.path.join(
                config["params"]["classify"]["kraken2"]["database"], "{db}"),
                              db = ["hash.k2d", "taxo.k2d", "opts.k2d"])
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
            quick = "--quick" \
                if config["params"]["classify"]["kraken2"]["quick"] \
                   else "",
            memory_mapping = "--memory-mapping" \
                if config["params"]["classify"]["kraken2"]["memory_mapping"] \
                   else "",
            use_names = "--use-names" \
                if config["params"]["classify"]["kraken2"]["use_names"] \
                   else "",
            use_mpa_style = "--use-mpa-style" \
                if config["params"]["classify"]["kraken2"]["use_mpa_style"] \
                   else "",
            report_zero_counts = "--report-zero-counts" \
                if config["params"]["classify"]["kraken2"]["report_zero_counts"] \
                   else "",
            confidence = config["params"]["classify"]["kraken2"]["confidence"],
            min_base_quality = config["params"]["classify"]["kraken2"]["min_base_quality"],
            min_hit_groups = config["params"]["classify"]["kraken2"]["min_hit_groups"],
            unclassified_out = "--unclassified-out %s" % \
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.unclassified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["classify"]["kraken2"]["unclassified_out"] \
                       else "",
            classified_out = "--classified-out %s" % \
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.classified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["classify"]["kraken2"]["classified_out"] \
                       else ""
        threads:
            config["params"]["classify"]["threads"]
        run:
            import os

            shell(
                '''
                kraken2 \
                {params.quick} \
                {params.memory_mapping} \
                {params.use_mpa_style} \
                {params.use_names} \
                {params.report_zero_counts} \
                --threads {threads} \
                --db {params.database} \
                --confidence {params.confidence} \
                --minimum-base-quality {params.min_base_quality} \
                --minimum-hit-groups {params.min_hit_groups} \
                {params.unclassified_out} \
                {params.classified_out} \
                --output %s \
                --report %s \
                --gzip-compressed \
                {params.paired} \
                {input.reads} \
                2> {log}
                ''' % (os.path.splitext(output.table)[0],
                       os.path.splitext(output.report)[0]))

            output_dir = os.path.dirname(output.report)
            with os.scandir(output_dir) as itr:
                for entry in itr:
                    if entry.is_file():
                        shell('''pigz -p {threads} %s''' \
                              % os.path.join(output_dir, entry.name))


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
