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


INDEX_SUFFIX = ["amb", "ann", "bwt", "pac", "sa"]
if config["params"]["alignment"]["aligner"] == "bwa-mem2":
    INDEX_SUFFIX = ["0123", "amb", "ann", "bwt.2bit.64", "pac"]
elif config["params"]["alignment"]["aligner"] == "bowtie2":
    INDEX_SUFFIX = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]


rule alignment_scaftigs_index:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "index/{{binning_group}}.{{assembly_group}}.{{assembler}}/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=INDEX_SUFFIX) if config["params"]["alignment"]["save_bam"] else \
        temp(expand(
            os.path.join(
                config["output"]["alignment"],
                "index/{{binning_group}}.{{assembly_group}}.{{assembler}}/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=INDEX_SUFFIX))
    log:
        os.path.join(config["output"]["alignment"], "logs/alignment_scaftigs_index/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(config["output"]["alignment"], "benchmark/alignment_scaftigs_index/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        aligner = config["params"]["alignment"]["aligner"],
        index = os.path.join(
            config["output"]["alignment"],
            "index/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"]
    conda:
        config["envs"]["align"]
    shell:
        '''
        if [ "{params.aligner}" == "bwa" ] || [ "{params.aligner}" == "bwa-mem2" ];
        then
            {params.aligner} index \
            {input} \
            -p {params.index} \
            >{log} 2>&1
        elif [ "{params.aligner}" == "bowtie2" ];
        then
            bowtie2-build \
            --threads {threads} \
            {input} \
            {params.index} \
            >{log} 2>&1
        fi
        '''


rule alignment_scaftigs_reads:
    input:
        reads = os.path.join(SAMPLESDIR, "reads/{sample}/{sample}.json"),
        index = expand(os.path.join(
            config["output"]["alignment"],
            "index/{{binning_group}}.{{assembly_group}}.{{assembler}}/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz.{suffix}"),
            suffix=INDEX_SUFFIX)
    output:
        stats = os.path.join(
            config["output"]["alignment"],
            "report/flagstat/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.flagstat"),
        bam = os.path.join(
            config["output"]["alignment"],
            "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam") \
            if config["params"]["alignment"]["save_bam"] else \
            temp(os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam")),
        bai = os.path.join(
            config["output"]["alignment"],
            "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam.bai") \
            if config["params"]["alignment"]["save_bam"] else \
            temp(os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam.bai"))
    log:
        os.path.join(
            config["output"]["alignment"],
            "logs/alignment_scaftigs_reads/{binning_group}.{assembly_group}.{assembler}/{sample}.log")
    benchmark:
        os.path.join(
            config["output"]["alignment"],
            "benchmark/alignment_scaftigs_reads/{binning_group}.{assembly_group}.{assembler}/{sample}.txt")
    params:
        aligner = config["params"]["alignment"]["aligner"],
        index = os.path.join(
            config["output"]["alignment"],
            "index/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    threads:
        config["params"]["alignment"]["threads"]
    conda:
        config["envs"]["align"]
    shell:
        '''
        OUTDIR=$(dirname {output.bam})
        OUTPE=$OUTDIR/pe
        OUTSE=$OUTDIR/se

        rm -rf $OUTDIR
        mkdir -p $OUTDIR

        R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
        R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
        RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

        STATSDIR=$(dirname {output.stats})
        rm -rf $STATSDIR
        mkdir -p $STATSDIR

        if [ "{params.aligner}" == "bwa" ] || [ "{params.aligner}" == "bwa-mem2" ];
        then
            if [ "$R1" != "" ];
            then
                mkdir -p $OUTPE
                STATSPE=$STATSDIR/{params.aligner}.pe.flagstat

                {params.aligner} mem \
                -t {threads} \
                {params.index} \
                $R1 $R2 \
                2> {log} | \
                tee >(samtools flagstat \
                -@4 - > $STATSPE) | \
                samtools sort \
                -m 3G \
                -@4 \
                -T $OUTPE/temp \
                -O BAM -o $OUTPE/sorted.bam -
            fi

            if [ "$RS" != "" ];
            then
                mkdir -p $OUTSE
                STATSSE=$STATSDIR/{params.aligner}.se.flagstat

                {params.aligner} mem \
                -t {threads} \
                {params.index} \
                $RS \
                2>> {log} | \
                tee >(samtools flagstat \
                -@4 - > $STATSSE) | \
                samtools sort \
                -m 3G \
                -@4 \
                -T $OUTSE/temp \
                -O BAM -o $OUTSE/sorted.bam -
            fi

            if [ -s $OUTPE/sorted.bam ] && [ -s $OUTSE/sorted.bam ];
            then
                samtools merge \
                -l 6 \
                -O BAM -o {output.bam} \
                $OUTPE/sorted.bam \
                $OUTSE/sorted.bam \
                2>>{log}

                rm -rf $OUTPE $OUTSE

            elif [ -s $OUTPE/sorted.bam ];
            then
                mv $OUTPE/sorted.bam {output.bam}
                mv $STATSPE {output.stats}
                rm -rf $OUTPE

            elif [ -s $OUTSE/sorted.bam ];
            then
                mv $OUTSE/sorted.bam {output.bam}
                mv $STATSSE {output.stats}
                rm -rf $OUTSE
            fi

            samtools flagstat \
            -@{threads} {output.bam} \
            > {output.stats} \
            2>>{log}

            samtools index \
            -@{threads} \
            {output.bam} {output.bai} \
            2>> {log}

        elif [ "{params.aligner}" == "bowtie2" ];
        then
            # see https://www.biostars.org/p/334422/

            READS=""
            if [ "$R1" != "" ] && [ "$RS" != "" ];
            then
                READS="-1 $R1 -2 $R2 -U $RS"
            elif [ "$R1" != "" ];
            then
                READS="-1 $R1 -2 $R2"
            elif [ "$RS" != "" ];
            then
                READS="-U $RS"
            fi

            bowtie2 \
            --threads {threads} \
            -x {params.index} \
            $READS \
            2> {log} | \
            tee >(samtools flagstat \
            -@4 - {output.stats}) | \
            samtools sort \
            -m 3G \
            -@4 \
            -T $OUTDIR/temp \
            -O BAM -o {output.bam} -

            rm -rf $OUTDIR/temp*

            samtools index \
            -@{threads} \
            {output.bam} {output.bai} \
            2>> {log}
        fi
        '''


rule alignment_scaftigs_reads_all:
    input:
        expand([
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.flagstat"),
            os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam"),
            os.path.join(
                config["output"]["alignment"],
                "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam.bai")],
            zip,
            binning_group=ALIGNMENT_GROUPS["binning_group"],
            assembly_group=ALIGNMENT_GROUPS["assembly_group"],
            assembler=ALIGNMENT_GROUPS["assembler"],
            sample=ALIGNMENT_GROUPS["sample_id"])


rule alignment_base_depth:
    input:
        os.path.join(
            config["output"]["alignment"],
            "bam/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.sorted.bam")
    output:
        os.path.join(
            config["output"]["alignment"],
            "depth/{binning_group}.{assembly_group}.{assembler}/{sample}/{sample}.align2scaftigs.depth.gz")
    log:
        os.path.join(
            config["output"]["alignment"],
            "logs/alignment_base_depth/{binning_group}.{assembly_group}.{assembler}/{sample}.log")
    benchmark:
        os.path.join(
            config["output"]["alignment"],
            "benchmark/alignment_base_depth/{binning_group}.{assembly_group}.{assembler}/{sample}.txt")
    threads:
        config["params"]["alignment"]["threads"]
    conda:
        config["envs"]["align"]
    shell:
        '''
        samtools depth {input} | \
        pigz -q -f -p {threads} -c \
        > {output} 2> {log}
        '''


if config["params"]["alignment"]["cal_base_depth"]:
    rule alignment_base_depth_all:
        input:
            expand(os.path.join(
                config["output"]["alignment"],
                "depth/{binning_group}.{assembly_group}.{assembler}/{sample}.align2scaftigs.depth.gz"),
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
                "report/flagstat/{binning_group}.{assembly_group}.{{assembler}}/{sample}/{sample}.align2scaftigs.flagstat"),
                zip,
                binning_group=ALIGNMENT_GROUP["binning_group"],
                assembly_group=ALIGNMENT_GROUP["assembly_group"],
                sample=ALIGNMENT_GROUP["sample_id"])
    output:
        flagstat = os.path.join(config["output"]["alignment"], "report/alignment_flagstat_{assembler}.tsv")
    run:
        input_list = [str(i) for i in input]
        metapi.flagstats_summary(input_list, 2, output=output.flagstat)


rule alignment_report_all:
    input:
        expand(os.path.join(config["output"]["alignment"], "report/alignment_flagstat_{assembler}.tsv"),
        assembler=ASSEMBLERS)


rule alignment_all:
    input:
        rules.alignment_base_depth_all.input,
        rules.alignment_report_all.input


localrules:
    alignment_scaftigs_reads_all,
    alignment_base_depth_all,
    alignment_report_all,
    alignment_all