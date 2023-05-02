rule raw_prepare_reads:
    input:
        lambda wildcards: metapi.get_raw_input_list(wildcards, SAMPLES, DATA_TYPE)
    output:
        os.path.join(config["output"]["raw"], "reads/{sample}/{sample}.json")
    log:
        os.path.join(config["output"]["raw"], "logs/raw_prepare_reads/{sample}.log")
    benchmark:
        os.path.join(config["output"]["raw"], "benchmark/raw_prepare_reads/{sample}.txt")
    params:
        sample_id = "{sample}",
        input_files = lambda wildcards: metapi.get_raw_input_dict(wildcards, SAMPLES, DATA_TYPE),
        headers = metapi.HEADERS,
        check_paired = config["params"]["raw"]["check_paired"]
    threads:
        config["params"]["raw"]["threads"]
    conda:
        config["envs"]["raw"]
    script:
        "../wrappers/preprocess_raw.py"


rule raw_prepare_reads_all:
    input:
        expand(os.path.join(config["output"]["raw"], "reads/{sample}/{sample}.json"),
        sample=SAMPLES_ID_LIST)


rule raw_fastqc:
    input:
        rules.raw_prepare_reads.output
    output:
        directory(os.path.join(config["output"]["raw"], "report/fastqc/{sample}.fastqc"))
    log:
        os.path.join(config["output"]["raw"], "logs/raw_fastqc/{sample}.log")
    benchmark:
        os.path.join(config["output"]["raw"], "benchmark/raw_fastqc/{sample}.txt")
    threads:
        config["params"]["raw"]["threads"]
    conda:
        config["envs"]["fastqc"]
    shell:
        '''
        mkdir -p {output}

        R1=$(jq -r -M '.PE_FORWARD' {input} | sed 's/^null$//g')
        R2=$(jq -r -M '.PE_REVERSE' {input} | sed 's/^null$//g')
        RS=$(jq -r -M '.SE' {input} | sed 's/^null$//g')

        fastqc \
        --outdir {output} \
        --threads {threads} \
        --format fastq \
        $R1 $R2 $RS \
        >{log} 2>&1
        '''


rule raw_fastqc_multiqc:
    input:
        expand(os.path.join(config["output"]["raw"], "report/fastqc/{sample}.fastqc"),
        sample=SAMPLES_ID_LIST)
    output:
        os.path.join(config["output"]["raw"], "report/multiqc/fastqc_multiqc_report.html")
    log:
        os.path.join(config["output"]["raw"], "logs/raw_fastqc_multiqc/raw_fastqc_multiqc.log")
    benchmark:
        os.path.join(config["output"]["raw"], "benchmark/raw_fastqc_multiqc/raw_fastqc_multiqc.txt")
    threads:
        1
    conda:
        config["envs"]["multiqc"]
    shell:
        '''
        OUTDIR=$(dirname {output})

        multiqc \
        --outdir $OUTDIR \
        --title fastqc \
        --module fastqc \
        {input} \
        >{log} 2>&1
        '''


rule raw_fastqc_all:
    input:
        expand([
            os.path.join(config["output"]["raw"], "report/fastqc/{sample}.fastqc"),
            os.path.join(config["output"]["raw"], "report/multiqc/fastqc_multiqc_report.html")],
            sample=SAMPLES_ID_LIST)


rule raw_report:
    input:
        rules.raw_prepare_reads.output
    output:
        os.path.join(config["output"]["raw"], "report/stats/{sample}_raw_stats.tsv")
    log:
        os.path.join(config["output"]["raw"], "logs/raw_report/{sample}.log")
    benchmark:
        os.path.join(config["output"]["raw"], "benchmark/raw_report/{sample}.txt")
    params:
        fq_encoding = config["params"]["fq_encoding"]
    threads:
        config["params"]["qcreport"]["seqkit"]["threads"]
    conda:
        config["envs"]["report"]
    shell:
        '''
        R1=$(jq -r -M '.PE_FORWARD' {input} | sed 's/^null$//g')
        R2=$(jq -r -M '.PE_REVERSE' {input} | sed 's/^null$//g')
        RS=$(jq -r -M '.SE' {input} | sed 's/^null$//g')
        RL=$(jq -r -M '.LONG' {input} | sed 's/^null$//g')

        seqkit stats \
        --all \
        --basename \
        --tabular \
        --fq-encoding {params.fq_encoding} \
        --out-file {output} \
        --threads {threads} \
        $R1 $R2 $RS $RL \
        >{log} 2>&1
        '''


rule raw_report_merge:
    input:
        expand(os.path.join(config["output"]["raw"], "report/stats/{sample}_raw_stats.tsv"),
        sample=SAMPLES_ID_LIST)
    output:
        os.path.join(config["output"]["qcreport"], "raw_stats.tsv")
    log:
        os.path.join(config["output"]["raw"], "logs/raw_report_merge/raw_report_merge.log")
    benchmark:
        os.path.join(config["output"]["raw"], "benchmark/raw_report_merge/raw_report_merge.txt")
    threads:
        config["params"]["qcreport"]["seqkit"]["threads"]
    shell:
        '''
        head -1 {input[0]} > {output}
        for report in {input}
        do
            tail -q -n +2 $report >> {output} 2>>{log}
        done
        '''


rule raw_report_all:
    input:
        os.path.join(config["output"]["qcreport"], "raw_stats.tsv")


rule raw_all:
    input:
        rules.raw_prepare_reads_all.input,
        rules.raw_fastqc_all.input,
        rules.raw_report_all.input


localrules:
    raw_prepare_reads_all,
    raw_fastqc_all,
    raw_report_all,
    raw_all
