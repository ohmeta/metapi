if config["params"]["trimming"]["sickle"]["do"]:
    rule trimming_sickle:
        input:
            rules.raw_prepare_reads.output
        output:
            os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["trimming"], "logs/trimming_sickle/{sample}.log")
        benchmark:
            os.path.join(config["output"]["trimming"], "benchmark/trimming_sickle/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["trimming"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["trimming"], "reads/{sample}/se/{sample}"),
            quality_type = config["params"]["trimming"]["sickle"]["quality_type"],
            quality_cutoff = config["params"]["trimming"]["sickle"]["quality_cutoff"],
            length_cutoff = config["params"]["trimming"]["sickle"]["length_cutoff"]
        priority:
            10
        threads:
            1
        conda:
            config["envs"]["trimming"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR

            R1=$(jq -r -M '.PE_FORWARD' {input} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""

            if [ $R1 != " " ];
            then
                FQ1={params.pe_prefix}.trimming.pe.1.fq.gz
                FQ2={params.pe_prefix}.trimming.pe.2.fq.gz
                FQPESI={params.pe_prefix}.trimming.pe.single.fq.gz

                mkdir -p $OUTPE

                sickle pe \
                --pe-file1 $R1 \
                --pe-file2 $R2 \
                --output-pe1 $FQ1 \
                --output-pe2 $FQ2 \
                --output-single $FQPESI \
                --qual-type {params.quality_type} \
                --qual-threshold {params.quality_cutoff} \
                --length-threshold {params.length_cutoff} \
                --gzip-output \
                >{log} 2>&1
            fi

            if [ $RS != " " ];
            then
                FQS={params.se_prefix}.trimming.se.fq.gz
                FQSESI={params.se_prefix}.trimming.se.single.fq.gz

                mkdir -p $OUTSE

                sickle se \
                --fastq-file $RS \
                --output-file $FQS \
                --output-single $FQSESI \
                --qual-type {params.quality_type} \
                --qual-threshold {params.quality_cutoff} \
                --length-threshold {params.length_cutoff} \
                --gzip-output \
                >>{log} 2>&1
            fi

            echo "{{\\"PE_FORWARD\\": \\"$FQ1\\", \\"PE_REVERSE\\": \\"$FQ2\\", \\"SE\\": \\"$FQS\\"}}" | \
            jq . > {output}
            '''


    rule trimming_sickle_all:
        input:
            expand(
                os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json"),
                sample=SAMPLES_ID_LIST)

else:
    rule trimming_sickle_all:
        input:


if config["params"]["trimming"]["fastp"]["do"]:
    rule trimming_fastp:
        input:
            rules.raw_prepare_reads.output
        output:
            os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["trimming"], "logs/trimming_fastp/{sample}.log")
        benchmark:
            os.path.join(config["output"]["trimming"], "benchmark/trimming_fastp/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["trimming"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["trimming"], "reads/{sample}/se/{sample}"),
            adapter_sequence = config["params"]["trimming"]["fastp"]["adapter_sequence"],
            adapter_sequence_r2 = config["params"]["trimming"]["fastp"]["adapter_sequence_r2"],
            adapter_operation = "--disable_adapter_trimming" \
            if config["params"]["trimming"]["fastp"]["disable_adapter_trimming"] else "",
            detect_adapter_pe = config["params"]["trimming"]["fastp"]["detect_adapter_for_pe"],
            detect_adapter_se = config["params"]["trimming"]["fastp"]["detect_adapter_for_se"],
            compression = config["params"]["trimming"]["fastp"]["compression"],
            cut_front_window_size = config["params"]["trimming"]["fastp"]["cut_front_window_size"],
            cut_front_mean_quality = config["params"]["trimming"]["fastp"]["cut_front_mean_quality"],
            cut_tail_window_size = config["params"]["trimming"]["fastp"]["cut_tail_window_size"],
            cut_tail_mean_quality = config["params"]["trimming"]["fastp"]["cut_tail_mean_quality"],
            cut_right_window_size = config["params"]["trimming"]["fastp"]["cut_right_window_size"],
            cut_right_mean_quality = config["params"]["trimming"]["fastp"]["cut_right_mean_quality"],
            length_required = config["params"]["trimming"]["fastp"]["length_required"],
            n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"],
            use_slide_window = "yes" if config["params"]["trimming"]["fastp"]["use_slide_window"] else "no",
            dedup = f'''--dedup --dup_calc_accuracy {config["params"]["trimming"]["fastp"]["dup_calc_accuracy"]}''' \
                if config["params"]["trimming"]["fastp"]["dedup"] else ""

        priority:
            10
        threads:
            config["params"]["trimming"]["threads"]
        conda:
            config["envs"]["trimming"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR

            R1=$(jq -r -M '.PE_FORWARD' {input} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""
            ADAPTER_OPERATION_PE=""
            ADAPTER_OPERATION_SE=""

            ADAPTER_OPERATION="{params.adapter_operation}"
            if [ $ADAPTER_OPERATION == "" ];
            then
                if [ "{params.detect_adapter_pe}" == "True" ]
                then
                    ADAPTER_OPERATION_PE="--detect_adapter_for_pe"
                else
                    ADAPTER_OPERATION_PE="--adapter_sequence {params.adapter_sequence} --adapter_sequence_r2 {params.adapter_sequence_r2}"
                fi

                if [ "{params.detect_adapter_se}" == "True" ]
                then
                    ADAPTER_OPERATION_SE="--detect_adapter_for_se"
                esle
                    ADAPTER_OPERATION_SE="--adapter_sequence {params.adapter_sequence}"
                fi
            fi

            if [ $R1 != "" ];
            then
                FQ1={params.pe_prefix}.trimming.pe.1.fq.gz
                FQ2={params.pe_prefix}.trimming.pe.2.fq.gz
                HTML={params.pe_prefix}.fastp.pe.html
                JSON={params.pe_prefix}.fastp.pe.json

                mkdir -p $OUTPE

                if [ "{params.use_slide_window}" == "yes" ];
                then
                    fastp \
                    --in1 $R1 \
                    --in2 $R2 \
                    --out1 $FQ1 \
                    --out2 $FQ2 \
                    --compression {params.compression} \
                    $ADAPTER_OPERATION_PE \
                    {params.dedup} \
                    --cut_front \
                    --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_right_window_size} \
                    --cut_right_mean_quality {params.cut_right_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html $HTML \
                    --json $JSON \
                    >{log} 2>&1
                else
                    fastp \
                    --in1 $R1 \
                    --in2 $R2 \
                    --out1 $FQ1 \
                    --out2 $FQ2 \
                    --compression {params.compression} \
                    $ADAPTER_OPERATION_PE \
                    {params.dedup} \
                    --cut_front \
                    --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html $HTML \
                    --json $JSON \
                    >{log} 2>&1
                fi
            fi

            if [ $RS != "" ];
            then
                FQS={params.se_prefix}.trimming.se.fq.gz
                HTML={params.se_prefix}.fastp.se.html
                JSON={params.se_prefix}.fastp.se.json

                mkdir -p $OUTSE

                if [ "{params.use_slide_window}" == "yes" ];
                then
                    fastp \
                    --in1 $RS \
                    --out1 $FQS \
                    --compression {params.compression} \
                    $ADAPTER_OPERATION_SE \
                    {params.dedup} \
                    --cut_front \
                    --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_right_window_size} \
                    --cut_right_mean_quality {params.cut_right_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html $HTML \
                    --json $JSON \
                    >>{log} 2>&1
                else
                    fastp \
                    --in1 $RS \
                    --out1 $FQS \
                    --compression {params.compression} \
                    $ADAPTER_OPERATION_SE \
                    {params.dedup} \
                    --cut_front \
                    --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html $HTML \
                    --json $JSON \
                    >>{log} 2>&1
                fi
            fi

            echo "{{\\"PE_FORWARD\\": \\"$FQ1\\", \\"PE_REVERSE\\": \\"$FQ2\\", \\"SE\\": \\"$FQS\\"}}" | \
            jq . > {output}
            '''


    rule trimming_fastp_all:
        input:
            expand(os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)

else:
    rule trimming_fastp_all:
        input:


if config["params"]["trimming"]["trimmomatic"]["do"]:
    rule trimming_trimmomatic:
        input:
            rules.raw_prepare_reads.output
        output:
            os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["trimming"], "logs/trimming_trimmomatic/{sample}.log")
        benchmark:
            os.path.join(config["output"]["trimming"], "benchmark/trimming_trimmomatic/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["trimming"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["trimming"], "reads/{sample}/se/{sample}"),
            trimmomatic_options = config["params"]["trimming"]["trimmomatic"]["trimmomatic_options"],
            phred = config["params"]["trimming"]["trimmomatic"]["phred"],
            save_unpaired = "true" if config["params"]["trimming"]["trimmomatic"]["save_unpaired"] else "false"
        priority:
            10
        threads:
            config["params"]["trimming"]["threads"]
        conda:
            config["envs"]["trimming"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR

            R1=$(jq -r -M '.PE_FORWARD' {input} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""

            if [ $R1 != "" ];
            then
                FQ1={params.pe_prefix}.trimming.pe.1.fq.gz
                FQ2={params.pe_prefix}.trimming.pe.2.fq.gz

                mkdir -p $OUTPE

                trimmomatic PE \
                {params.phred} \
                -threads {threads} \
                -summary {params.pe_prefix}.summary.txt \
                $R1 $R2 \
                $FQ1 $FQ1.unpaired.gz \
                $FQ2 $FQ2.unpaired.gz \
                {params.trimmomatic_options} \
                >{log} 2>&1

                if [ "{params.save_unpaired}" == "false" ];
                then
                    rm -rf $FQ1.unpaired.gz
                    rm -rf $FQ2.unpaired.gz
                fi
            fi

            if [ $RS != "" ];
            then
                FQS={params.se_prefix}.trimming.se.fq.gz

                mkdir -p $OUTSE

                trimmomatic SE \
                {params.phred} \
                -threads {threads} \
                -summary {params.se_prefix}.summary.txt \
                $RS \
                $FQS \
                {params.trimmomatic_options} \
                >>{log} 2>&1
            fi

            echo "{{\\"PE_FORWARD\\": \\"$FQ1\\", \\"PE_REVERSE\\": \\"$FQ2\\", \\"SE\\": \\"$FQS\\"}}" | \
            jq . > {output}
            '''


    rule trimming_trimmomatic_all:
        input:
            expand(os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)

else:
    rule trimming_trimmomatic_all:
        input:


if TRIMMING_DO and config["params"]["qcreport"]["do"]:
    rule trimming_report:
        input:
            os.path.join(config["output"]["trimming"], "reads/{sample}/{sample}.json")
        output:
            os.path.join(config["output"]["trimming"], "report/stats/{sample}_trimming_stats.tsv")
        log:
            os.path.join(config["output"]["trimming"], "logs/trimming_report/{sample}.log")
        benchmark:
            os.path.join(config["output"]["trimming"], "benchmark/trimming_report/{sample}.txt")
        params:
            fq_encoding = config["params"]["fq_encoding"]
        priority:
            10
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        conda:
            config["envs"]["report"]
        shell:
            '''
            R1=$(jq -r -M '.PE_FORWARD' {input} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input} | sed 's/^null$//g')

            seqkit stats \
            --all \
            --basename \
            --tabular \
            --fq-encoding {params.fq_encoding} \
            --out-file {output} \
            --threads {threads} \
            $R1 $R2 $RS \
            >{log} 2>&1
            '''


    rule trimming_report_merge:
        input:
            expand(os.path.join(config["output"]["trimming"], "report/stats/{sample}_trimming_stats.tsv"),
            sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "trimming_stats.tsv")
        log:
            os.path.join(config["output"]["raw"], "logs/trimming_report_merge/trimming_report_merge.log")
        benchmark:
            os.path.join(config["output"]["raw"], "benchmark/trimming_report_merge/trimming_report_merge.txt")
        priority:
            10
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


    rule trimming_report_all:
        input:
            rules.trimming_report_merge.output

else:
    rule trimming_report_all:
        input:


rule trimming_all:
    input:
        rules.trimming_sickle_all.input,
        rules.trimming_fastp_all.input,
        rules.trimming_trimmomatic_all.input,
        rules.trimming_report_all.input


localrules:
    trimming_fastp_all,
    trimming_sickle_all,
    trimming_trimmomatic_all,
    trimming_report_all,
    trimming_all
