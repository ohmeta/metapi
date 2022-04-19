ALIGNMENT_GROUP = SAMPLES.reset_index().loc[:, ["sample_id", "assembly_group"]].drop_duplicates()

alignment_df_list = []
for assembler in ASSEMBLERS:
    alignment_df = ALIGNMENT_GROUP.copy()
    alignment_df["assembler"] = assembler
    alignment_df_list.append(alignment_df)
ALIGNMENT_GROUPS = pd.concat(alignment_df_list, axis=0)


def alignment_input_with_short_reads(wildcards):
    if RMHOST_DO:
        return get_reads(wildcards, "rmhost", False, False)
    elif TRIMMING_DO:
        return get_reads(wildcards, "trimming", False, False)
    else:
        return get_reads(wildcards, "raw", False, False)


rule alignment_scaftigs_index:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        temp(expand(
            os.path.join(
                config["output"]["alignment"],
                "index/{{assembly_group}}.{{assembler}}.out/{{assembly_group}}.scaftigs.fa.gz.{suffix}"),
            suffix=BWA_INDEX_SUFFIX))
    log:
        os.path.join(
            config["output"]["alignment"],
            "logs/index/{assembly_group}.{assembler}.scaftigs.index.log")
    params:
        bwa = "bwa-mem2" if config["params"]["alignment"]["algorithms"] == "mem2" else "bwa",
        output_prefix = os.path.join(
            config["output"]["alignment"],
            "index/{assembly_group}.{assembler}.out/{assembly_group}.scaftigs.fa.gz")
    shell:
        '''
        {params.bwa} index {input.scaftigs} -p {params.output_prefix} 2> {log}
        '''


rule alignment_reads_scaftigs:
    input:
        reads = alignment_input_with_short_reads,
        index = expand(os.path.join(
            config["output"]["alignment"],
            "index/{{assembly_group}}.{{assembler}}.out/{{assembly_group}}.scaftigs.fa.gz.{suffix}"),
            suffix=BWA_INDEX_SUFFIX)
    output:
        flagstat = os.path.join(
            config["output"]["alignment"],
            "report/flagstat/{assembly_group}.{assembler}/{sample}.align2scaftigs.flagstat"),
        bam = temp(os.path.join(
            config["output"]["alignment"],
            "bam/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam")),
        bai = temp(os.path.join(
            config["output"]["alignment"],
            "bam/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam.bai"))
    log:
        os.path.join(config["output"]["alignment"],
                     "logs/alignment/{assembly_group}.{assembler}/{sample}.align2scaftigs.log")
    benchmark:
        os.path.join(config["output"]["alignment"],
                     "benchmark/alignment/{assembly_group}.{assembler}/{sample}.align2scaftigs.benchmark.txt")
    params:
        bwa = "bwa-mem2" if config["params"]["alignment"]["algorithms"] == "mem2" else "bwa",
        index_prefix = os.path.join(
            config["output"]["alignment"],
            "index/{assembly_group}.{assembler}.out/{assembly_group}.scaftigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        rm -rf {output.bam}*

        {params.bwa} mem \
        -t {threads} \
        {params.index_prefix} \
        {input.reads} 2> {log} |
        tee >(samtools flagstat \
              -@{threads} - \
              > {output.flagstat}) | \
        samtools sort \
        -m 3G \
        -@{threads} \
        -T {output.bam} \
        -O BAM -o {output.bam} -

        samtools index -@{threads} {output.bam} {output.bai} 2>> {log}
        '''


rule alignment_reads_scaftigs_all:
    input:
        expand([
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat/{assembly_group}.{assembler}/{sample}.align2scaftigs.flagstat"),
            os.path.join(
                config["output"]["alignment"],
                "bam/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam"),
            os.path.join(
                config["output"]["alignment"],
                "bam/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam.bai")],
            zip,
            assembly_group=ALIGNMENT_GROUPS["assembly_group"],
            assembler=ALIGNMENT_GROUPS["assembler"],
            sample=ALIGNMENT_GROUPS["sample_id"])


if config["params"]["alignment"]["cal_base_depth"]:
    rule alignment_base_depth:
        input:
            os.path.join(
                config["output"]["alignment"],
                "bam/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam")
        output:
            os.path.join(
                config["output"]["alignment"],
                "depth/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.depth.gz")
        shell:
            '''
            samtools depth {input} | gzip -c > {output}
            '''


    rule alignment_base_depth_all:
        input:
            expand(os.path.join(
                config["output"]["alignment"],
                "depth/{assembly_group}.{assembler}.out/{sample}.align2scaftigs.depth.gz"),
                zip,
                assembly_group=ALIGNMENT_GROUPS["assembly_group"],
                assembler=ALIGNMENT_GROUPS["assembler"],
                sample=ALIGNMENT_GROUPS["sample_id"])

else:
    rule alignment_base_depth_all:
        input:


rule alignment_report:
    input:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat/{assembly_group}.{{assembler}}/{sample}.align2scaftigs.flagstat"),
                zip,
                assembly_group=ALIGNMENT_GROUP["assembly_group"],
                sample=ALIGNMENT_GROUP["sample_id"])
    output:
        flagstat = os.path.join(config["output"]["alignment"],
                                "report/alignment_flagstat_{assembler}.tsv")
    run:
        input_list = [str(i) for i in input]
        metapi.flagstats_summary(input_list, 2, output=output.flagstat)


rule alignment_report_all:
    input:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "report/alignment_flagstat_{assembler}.tsv"),
            assembler=ASSEMBLERS)


rule alignment_all:
    input:
        #rules.alignment_reads_scaftigs_all,
        rules.alignment_base_depth_all.input,
        rules.alignment_report_all.input


localrules:
    alignment_base_depth_all,
    alignment_reads_scaftigs_all,
    alignment_report,
    alignment_report_all,
    alignment_all