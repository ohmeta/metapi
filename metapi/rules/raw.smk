def raw_input(wildcards):
    fqs = {}
    headers = HEADERS[DT]
    for k, v in headers.items():
        if v in SAMPLES.columns:
            reads = metapi.get_reads(SAMPLES, wildcards, v)
            if len(reads) > 0:
                fqs[v] = reads
    return fqs


rule prepare_reads:
    input:
        unpack(raw_input)
    output:
        os.path.join(config["output"]["raw"], "reads/{sample}/done")
    params:
        sample_id = "{sample}",
        input_files = lambda wildcards: raw_input(wildcards),
        headers = HEADERS,
        check_paired = config["params"]["raw"]["check_paired"]
    threads:
        config["params"]["raw"]["threads"]
    conda:
        config["envs"]["raw"]
    log:
        os.path.join(config["output"]["raw"], "logs/{sample}_prepare.log")
    script:
        "../wrappers/preprocess_raw.py"


rule prepare_reads_all:
    input:
        expand(os.path.join(config["output"]["raw"], "reads/{sample}/done"),
        sample=SAMPLES_ID_LIST)


rule raw_fastqc:
    input:
        rules.prepare_reads.output
    output:
        directory(os.path.join(config["output"]["raw"], "fastqc/{sample}.fastqc.out"))
    conda:
        config["envs"]["fastqc"]
    threads:
        config["params"]["raw"]["threads"]
    log:
        os.path.join(config["output"]["raw"], "logs/{sample}.fastqc.log")
    shell:
        '''
        mkdir -p {output}

        fastqc \
        --outdir {output} \
        --threads {threads} \
        --format fastq \
        {input} \
        2> {log}
        '''


rule raw_fastqc_multiqc:
    input:
        expand(os.path.join(config["output"]["raw"], "fastqc/{sample}.fastqc.out"),
        sample=SAMPLES_ID_LIST)
    output:
        html = os.path.join(config["output"]["raw"], "report/fastqc_multiqc_report.html")
    conda:
        config["envs"]["multiqc"]
    params:
        outdir = os.path.join(config["output"]["raw"], "report")
    log:
        os.path.join(config["output"]["raw"], "logs/multiqc_fastqc.log")
    threads:
        1
    shell:
        '''
        multiqc \
        --outdir {params.outdir} \
        --title fastqc \
        --module fastqc {input} \
        2> {log}
        '''


rule raw_fastqc_all:
    input:
        expand([
            os.path.join(config["output"]["raw"], "fastqc/{sample}.fastqc.out"),
            os.path.join(config["output"]["raw"], "report/fastqc_multiqc_report.html")],
            sample=SAMPLES_ID_LIST)


rule raw_report:
    input:
        rules.prepare_reads.output
    output:
        temp(os.path.join(config["output"]["raw"], "report/stats/{sample}_raw_stats.tsv.raw"))
    conda:
        config["envs"]["report"]
    params:
        fq_encoding = config["params"]["fq_encoding"]
    log:
        os.path.join(config["output"]["raw"], "logs/{sample}.seqkit.log")
    threads:
        config["params"]["qcreport"]["seqkit"]["threads"]
    shell:
        '''
        seqkit stats \
        --all \
        --basename \
        --tabular \
        --fq-encoding {params.fq_encoding} \
        --out-file {output} \
        --threads {threads} \
        {input} 2> {log}
        '''


rule raw_report_refine:
    input:
        os.path.join(config["output"]["raw"], "report/stats/{sample}_raw_stats.tsv.raw")
    output:
        os.path.join(config["output"]["raw"], "report/stats/{sample}_raw_stats.tsv.gz")
    params:
        sample_id = "{sample}"
    threads:
        1
    run:
        if IS_PE:
            metapi.change(str(input), str(output), params.sample_id, "raw",
                            "pe", ["fq1", "fq2"])
        else:
            metapi.change(str(input), str(output), params.sample_id, "raw",
                            "se", ["fq1"])


rule raw_report_merge:
    input:
        expand(os.path.join(config["output"]["raw"], "report/stats/{sample}_raw_stats.tsv.gz"),
        sample=SAMPLES_ID_LIST)
    output:
        os.path.join(config["output"]["qcreport"], "raw_stats.tsv.gz")
    threads:
        config["params"]["qcreport"]["seqkit"]["threads"]
    run:
        metapi.merge(input, metapi.parse, threads, output=output[0])


rule raw_report_all:
    input:
        os.path.join(config["output"]["qcreport"], "raw_stats.tsv.gz")


rule raw_all:
    input:
        rules.prepare_reads_all.input,
        rules.raw_fastqc_all.input,
        rules.raw_report_all.input


localrules:
    prepare_reads_all,
    raw_fastqc_all,
    raw_report_all,
    raw_all
