rule classify_short_reads_kraken2:
    input:
        assembly_input
    output:
        table = os.path.join(
            config["output"]["classify"],
            "short_reads/{sample}.kraken2.out/{sample}.kraken2.table"),
        report = os.path.join(
            config["output"]["classify"],
            "short_reads/{sample}.kraken2.out/{sample}.kraken2.report")
    log:
        os.path.join(config["output"]["classify"],
                     "logs/{sample}.kraken2.log")
    params:
        paired = "--paired" if IS_PE else "",
        database = config["params"]["classify"]["kraken2"]["database"],
        report_zero_counts = \
            config["params"]["classify"]["kraken2"]["report_zero_counts"]
    threads:
        config["params"]["classify"]["kraken2"]["threads"]
    shell:
        '''
        kraken2 \
        --use-names \
        --threads {threads} \
        --db {params.database} \
        --output {output.table} \
        --report {output.report} \
        {params.report_zero_counts} \
        {params.paired} \
        --gzip-compressed \
        {input} \
        2> {log}
        '''


rule classify_short_reads_kraken2_all:
    input:
        expand([
            os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.table"),
            os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.report")],
               sample=SAMPLES.index.unique())


rule classify_hmq_bins_gtdbtk:
    input:
        os.path.join(
            config["output"]["checkm"],
            "hmq_bins/{assembler}.{binner}.links")
    output:
        directory(os.path.join(
            config["output"]["classify"],
            "hmq_bins/{assembler}.{binner}.gtdbtk.out"))
    log:
        os.path.join(config["output"]["classify"],
                     "logs/{assembler}.{binner}.gtdbtk.log")
    params:
        bin_suffix = "fa"
    threads:
        config["params"]["classify"]["gtdbtk"]["threads"]
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
                "hmq_bins/{assembler}.{binner}.gtdbtk.out"),
            assembler=ASSEMBLERS,
            binner=BINNERS)


rule classify_all:
    input:
        rules.classify_short_reads_kraken2_all.input,
        rules.classify_hmq_bins_gtdbtk_all.input
