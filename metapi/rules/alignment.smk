rule build_asmfa_index:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    output:
        expand("{assembly}/{{sample}}.megahit_out/{{sample}}.contigs.fa.gz.{suffix}",
               assembly=config["results"]["assembly"],
               suffix=["amb", "ann", "bwt", "pac", "sa"])
    log:
        os.path.join(config["logs"]["alignment"], "{sample}_bwa_index.log")
    shell:
        "bwa index {input} 2> {log}"


rule align_reads_to_asmfa:
    input:
        reads = clean_reads,
        index = expand("{assembly}/{{sample}}.megahit_out/{{sample}}.contigs.fa.gz.{suffix}",
                       assembly=config["results"]["assembly"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["alignment"], "{sample}.flagstat"),
        bam = os.path.join(config["results"]["alignment"], "{sample}.sorted.bam")
    log:
        os.path.join(config["logs"]["alignment"], "{sample}_bwa_mem.log")
    params:
        prefix = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"],
    shell:
        "bwa mem -t {threads} {params.prefix} {input.reads} | "
        "samtools view -@{threads} -hbS - | "
        "tee >(samtools flagstat -@{threads} - > {output.flagstat}) | "
        "samtools sort -@{threads} -o {output.bam} - 2> {log}"

rule index_bam:
    input:
        os.path.join(config["results"]["alignment"], "{sample}.sorted.bam")
    output:
        os.path.join(config["results"]["alignment"], "{sample}.sorted.bam.bai")
    log:
        os.path.join(config["logs"]["alignment"], "{sample}_bam_index.log")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        samtools index -@{threads} {input} {output}
        '''
