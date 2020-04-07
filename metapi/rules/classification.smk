rule kraken2:
    input:
        clean_reads
    output:
        out = os.path.join(config["results"]["classification"]["kraken2"], "{sample}.kraken2.output"),
        report = os.path.join(config["results"]["classification"]["kraken2"], "{sample}.kraken2.report")
    params:
        database = config["params"]["classification"]["kraken2"]["database"]
    threads:
        config["params"]["classification"]["kraken2"]["threads"]
    log:
        os.path.join(config["logs"]["classification"]["kraken2"], "{sample}.kraken2.log")
    shell:
        '''
        kraken2 \
        --use-names \
        --threads {threads} \
        --db {params.database} \
        --output {output.out} \
        --report {output.report} \
        --report-zero-counts \
        --paired \
        --gzip-compressed \
        {input[0]} {input[1]} \
        2> {log}
        '''


rule classification_hmq_bins_by_gtdbtk:
    input:
        bins_hmq = directory(os.path.join(config["results"]["checkm"]["base_dir"],
                                          "bins.{assembler}.{binner}_out.hmq"))
    output:
        outdir = directory(config["results"]["classification"]["gtdbtk"]["out_dir"])
    params:
        extension = config["params"]["classification"]["gtdbtk"]["extension"],
        scratch_dir = "--scratch_dir %s" % config["results"]["classification"]["gtdbtk"]["scratch_dir"] if config["params"]["classification"]["gtdbtk"]["reduce_memory"] else ""
    threads:
        config["params"]["classification"]["gtdbtk"]["threads"]
    log:
        os.path.join(config["logs"]["classification"]["gtdbtk"], "bins.hmq.gtdbtk.log")
    shell:
        '''
        rm -rf {output.outdir}
        rm -rf {params.scratch_dir}

        gtdbtk classify_wf \
        --genome_dir {input.bins_hmq}/ \
        --out_dir {output.outdir} \
        --extension {params.extension} \
        --cpus {threads} \
        {params.scratch_dir} >{log}
        '''
