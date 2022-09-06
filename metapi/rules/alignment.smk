ALIGNMENT_GROUP = SAMPLES.reset_index().loc[:, ["sample_id", "assembly_group", "binning_group"]].drop_duplicates()

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
            "scaftigs/{binning_group}.{assembly_group}.{assembler}.out/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "index/{{binning_group}}.{{assembly_group}}.{{assembler}}.out/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=BWA_INDEX_SUFFIX) if config["params"]["alignment"]["save_bam"] else \
        temp(expand(
            os.path.join(
                config["output"]["alignment"],
                "index/{{binning_group}}.{{assembly_group}}.{{assembler}}.out/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=BWA_INDEX_SUFFIX))
    conda:
        config["envs"]["align"]
    log:
        os.path.join(
            config["output"]["alignment"],
            "logs/index/{binning_group}.{assembly_group}.{assembler}.scaftigs.index.log")
    params:
        bwa = "bwa-mem2" if config["params"]["alignment"]["algorithms"] == "mem2" else "bwa",
        output_prefix = os.path.join(
            config["output"]["alignment"],
            "index/{binning_group}.{assembly_group}.{assembler}.out/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    shell:
        '''
        {params.bwa} index {input.scaftigs} -p {params.output_prefix} 2> {log}
        '''


rule alignment_reads_scaftigs:
    input:
        reads = alignment_input_with_short_reads,
        index = expand(os.path.join(
            config["output"]["alignment"],
            "index/{{binning_group}}.{{assembly_group}}.{{assembler}}.out/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=BWA_INDEX_SUFFIX)
    output:
        flagstat = os.path.join(
            config["output"]["alignment"],
            "report/flagstat/{binning_group}.{assembly_group}.{assembler}/{sample}.align2scaftigs.flagstat"),
        bam = os.path.join(
            config["output"]["alignment"],
            "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam") \
            if config["params"]["alignment"]["save_bam"] else \
            temp(os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam")),
        bai = os.path.join(
            config["output"]["alignment"],
            "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam.bai") \
            if config["params"]["alignment"]["save_bam"] else \
            temp(os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam.bai"))
    conda:
        config["envs"]["align"]
    log:
        os.path.join(config["output"]["alignment"],
                     "logs/alignment/{binning_group}.{assembly_group}.{assembler}/{sample}.align2scaftigs.log")
    benchmark:
        os.path.join(config["output"]["alignment"],
                     "benchmark/alignment/{binning_group}.{assembly_group}.{assembler}/{sample}.align2scaftigs.benchmark.txt")
    params:
        bwa = "bwa-mem2" if config["params"]["alignment"]["algorithms"] == "mem2" else "bwa",
        index_prefix = os.path.join(
            config["output"]["alignment"],
            "index/{binning_group}.{assembly_group}.{assembler}.out/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
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
                "report/flagstat/{binning_group}.{assembly_group}.{assembler}/{sample}.align2scaftigs.flagstat"),
            os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam"),
            os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam.bai")],
            zip,
            binning_group=ALIGNMENT_GROUPS["binning_group"],
            assembly_group=ALIGNMENT_GROUPS["assembly_group"],
            assembler=ALIGNMENT_GROUPS["assembler"],
            sample=ALIGNMENT_GROUPS["sample_id"])


if config["params"]["alignment"]["cal_base_depth"]:
    rule alignment_base_depth:
        input:
            os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.sorted.bam")
        output:
            os.path.join(
                config["output"]["alignment"],
                "depth/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.depth.gz")
        conda:
            config["envs"]["align"]
        shell:
            '''
            samtools depth {input} | gzip -c > {output}
            '''


    rule alignment_base_depth_all:
        input:
            expand(os.path.join(
                config["output"]["alignment"],
                "depth/{binning_group}.{assembly_group}.{assembler}.out/{sample}.align2scaftigs.depth.gz"),
                zip,
                binning_group=ALIGNMENT_GROUPS["binning_group"],
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
                "report/flagstat/{binning_group}.{assembly_group}.{{assembler}}/{sample}.align2scaftigs.flagstat"),
                zip,
                binning_group=ALIGNMENT_GROUP["binning_group"],
                assembly_group=ALIGNMENT_GROUP["assembly_group"],
                sample=ALIGNMENT_GROUP["sample_id"])
    output:
        flagstat = os.path.join(config["output"]["alignment"],
                                "report/alignment_flagstat_{assembler}_bwa.tsv")
    run:
        input_list = [str(i) for i in input]
        metapi.flagstats_summary(input_list, 2, output=output.flagstat)


rule alignment_report_all:
    input:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "report/alignment_flagstat_{assembler}_bwa.tsv"),
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