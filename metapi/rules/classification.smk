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