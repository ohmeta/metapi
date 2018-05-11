def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["1.fq.gz", "2.fq.gz"]].dropna()

rule trimming_pe:
    input:
        fastq1=get_fastq[1]
        fastq2=
    output:
        fastq1="data/01.trimmed/{sample}_{unit}.trimmed.1.fq.gz",
        fastq2="data/01.trimmed/{sample}_{unit}.trimmed.2.fq.gz"
    log:
    params:
    wrapper:
        0.20/bio/sickle


rule trimming_se:
    input:
        fastq1=
        fastq2=
    output:
        fastq="data/01.trimmed/{sample}_{unit}.trimmed.fq.gz"
    log:
    params:
    wrapper:
