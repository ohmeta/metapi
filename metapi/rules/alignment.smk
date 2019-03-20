rule build_index_for_scaftigs:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        expand("{assembly}/{{sample}}.{{assembler}}_out/{{sample}}.{{assembler}}.scaftigs.fa.gz.{suffix}",
               assembly=config["results"]["assembly"],
               suffix=["amb", "ann", "bwt", "pac", "sa"])
    params:
        prefix = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.{assembler}_build_index_for_scaftigs.log")
    shell:
        '''
        bwa index {input} -p {params.prefix} 2> {log}
        '''


rule align_reads_to_scaftigs:
    input:
        reads = clean_reads,
        index = expand("{assembly}/{{sample}}.{{assembler}}_out/{{sample}}.{{assembler}}.scaftigs.fa.gz.{suffix}",
                       assembly=config["results"]["assembly"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.flagstat"),
        bam = temp(os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam"))
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.{assembler}_align_reads_to_scaftigs.log")
    params:
        prefix = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        bwa mem -t {threads} {params.prefix} {input.reads} |
        tee >(samtools flagstat -@ {threads} - > {output.flagstat}) |
        samtools sort -@{threads} -O BAM -o {output.bam} - 2> {log}
        '''


rule build_index_for_bam:
    input:
        os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam")
    output:
        temp(os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam.bai"))
    log:
        os.path.join(config["logs"]["alignment"], "{sample}.{assembler}_build_index_for_bam.log")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        samtools index -@{threads} {input} {output}
        '''


rule summary_scaftigs_flagstat:
    input:
        expand(os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.flagstat"),
               sample=_samples.index,
               assembler=config["params"]["assembler"])
    output:
        os.path.join(config["results"]["alignment"], "scaftigs_flagstat_summary.tsv")
    run:
        input_list = []
        for i in input:
            input_list.append(str(i))
        output_str = str(output)
        flagstats_summary(input_list, output_str, 2)
