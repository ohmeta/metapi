rule bwa_index_scaftigs:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        expand("{assembly}/{{sample}}.{{assembler}}_out/{{sample}}.{{assembler}}.scaftigs.fa.gz.{suffix}",
               assembly=config["results"]["assembly"],
               suffix=["amb", "ann", "bwt", "pac", "sa"])
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.{assembler}_bwa_index_scaftigs.log")
    shell:
        "bwa index {input} 2> {log}"


rule bwa_mem_scaftigs:
    input:
        reads = clean_reads,
        index = expand("{assembly}/{{sample}}.{{assembler}}_out/{{sample}}.{{assembler}}.scaftigs.fa.gz.{suffix}",
                       assembly=config["results"]["assembly"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.flagstat"),
        bam = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam")
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.{assembler}_bwa_mem_scaftigs.log")
    params:
        prefix = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"],
    shell:
        "bwa mem -t {threads} {params.prefix} {input.reads} | "
        "tee >(samtools flagstat -@ {threads} - > {output.flagstat}) | "
        "samtools sort -@ {threads} -o {output.bam} - 2> {log}"

rule bwa_index_bam:
    input:
        os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam")
    output:
        os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam.bai")
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.{assembler}_bam_index_bam.log")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        samtools index -@ {threads} {input} {output}
        '''
