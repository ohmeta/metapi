rule alignment_scaftigs_index:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "index/{{sample}}.{{assembler}}.out/{{sample}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=["amb", "ann", "bwt", "pac", "sa"])
    log:
        os.path.join(
            config["output"]["alignment"],
            "logs/index/{sample}.{assembler}.scaftigs.index.log")
    params:
        output_prefix = os.path.join(
            config["output"]["alignment"],
            "index/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
    shell:
        '''
        bwa index {input.scaftigs} -p {params.output_prefix} 2> {log}
        '''


rule alignment_reads_scaftigs:
    input:
        reads = assembly_input,
        index = expand(os.path.join(
            config["output"]["alignment"],
            "index/{{sample}}.{{assembler}}.out/{{sample}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = protected(
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat/{sample}.{assembler}.align2scaftigs.flagstat")),
        bam = temp(
            os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam"))
    log:
        os.path.join(config["output"]["alignment"],
                     "logs/alignment/{sample}.{assembler}.align.reads2scaftigs.log")
    params:
        index_prefix = os.path.join(
            config["output"]["alignment"],
            "index/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        bwa mem \
        -t {threads} \
        {params.index_prefix} \
        {input.reads} 2> {log} |
        tee >(samtools flagstat \
              -@{threads} - \
              > {output.flagstat}) | \
        samtools sort \
        -@{threads} \
        -T {output.bam} \
        -O BAM -o {output.bam} -
        '''


rule alignment_bam_index:
    input:
        os.path.join(
            config["output"]["alignment"],
            "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam")
    output:
        temp(
            os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam.bai"))
    log:
        os.path.join(config["output"]["alignment"],
                     "logs/{sample}.{assembler}.bam.index.log")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        samtools index -@{threads} {input} {output} 2> {log}
        '''

if config["params"]["alignment"]["cal_base_depth"]:
    rule alignment_base_depth:
        input:
            os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam")
        output:
            protected(os.path.join(
                config["output"]["alignment"],
                "depth/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.depth.gz"))
        shell:
            '''
            samtools depth {input} | gzip -c > {output}
            '''


    rule alignment_base_depth_all:
        input:
            expand(os.path.join(
                config["output"]["alignment"],
                "depth/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.depth.gz"),
                   assembler=ASSEMBLERS,
                   sample=SAMPLES.index.unique())

else:
    rule alignment_base_depth_all:
        input:


rule alignment_report:
    input:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat/{sample}.{{assembler}}.align2scaftigs.flagstat"),
            sample=SAMPLES.index.unique())
    output:
        os.path.join(config["output"]["alignment"],
                     "report/alignment_flagstat_{assembler}.tsv")
    run:
        input_list = [str(i) for i in input]
        output_str = str(output)
        metapi.flagstats_summary(input_list, output_str, 2)


rule alignment_all:
    input:
        expand(os.path.join(
            config["output"]["alignment"],
            "report/alignment_flagstat_{assembler}.tsv"),
               assembler=ASSEMBLERS,
               sample=SAMPLES.index.unique()),
        rules.alignment_base_depth_all.input,

        rules.assembly_all.input
