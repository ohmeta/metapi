rule burst_reads:
    input:
        reads = clean_reads
    output:
        align = os.path.join(config["results"]["burst"], "{sample}.reads.burst.b6"),
        tmp = temp(directory(os.path.join(config["results"]["burst"], "{sample}.tmp")))
    params:
        references = config["params"]["burst"]["references"],
        accelerator = config["params"]["burst"]["accelerator"],
        identity = config["params"]["burst"]["identity"],
        mode = config["params"]["burst"]["mode"],
        prefix = "{sample}",
        outdir = directory(config["results"]["burst"])
    log:
        os.path.join(config["logs"]["burst"], "{sample}.reads.burst.log")
    threads:
        config["params"]["burst"]["threads"]
    shell:
        '''
        mkdir -p {output.tmp}
        gzip -d {input.reads[0]} -c > {output.tmp}/{params.prefix}.trimmed.1.fq
        gzip -d {input.reads[1]} -c > {output.tmp}/{params.prefix}.trimmed.2.fq
        fq2fa --merge {output.tmp}/{params.prefix}.trimmed.1.fq {output.tmp}/{params.prefix}.trimmed.2.fq {output.tmp}/{params.prefix}.trimmed.pe.fa
        rm -rf {output.tmp}/{params.prefix}.trimmed.1.fq {output.tmp}/{params.prefix}.trimmed.2.fq
        burst12 --references {params.references} --accelerator {params.accelerator} \
        --queries {output.tmp}/{params.prefix}.trimmed.pe.fa \
        -o {output.align} \
        --threads {threads} \
        --mode {params.mode} >{log} 2>&1 &
        '''

# rule burst_contigs:

# rule burst_bins:
