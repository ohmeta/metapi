def rmhost_input(wildcards):
    if TRIMMING_DO:
        return os.path.join(config["output"]["trimming"], f"reads/{wildcards.sample}/{wildcards.sample}.json")
    else:
        return os.path.join(config["output"]["raw"], f"reads/{wildcards.sample}/{wildcards.sample}.json")


BWA_INDEX_SUFFIX = ["0123", "amb", "ann", "bwt.2bit.64", "pac"]
if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2":
    BWA_INDEX_SUFFIX = ["amb", "ann", "bwt", "pac", "sa"]


if config["params"]["rmhost"]["bwa"]["do"]:
    rule rmhost_bwa_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}", prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
            suffix=BWA_INDEX_SUFFIX)
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_bwa_index/rmhost_bwa_index.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_bwa_index/rmhost_bwa_index.txt")
        params:
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bwa = "bwa-mem2" if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" else "bwa"
        priority:
            20
        conda:
            config["envs"]["align"]
        shell:
            '''
            {params.bwa} index -p {params.index_prefix} {input} >{log} 2>&1
            '''


    rule rmhost_bwa:
        input:
            reads = lambda wildcards: rmhost_input(wildcards),
            index = expand(
                "{prefix}.{suffix}",
                prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                suffix=BWA_INDEX_SUFFIX)
        output:
            os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_bwa/{sample}.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_bwa/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/se/{sample}"),
            report_dir = os.path.join(config["output"]["rmhost"], "report/flagstat/{sample}"),
            compression = config["params"]["rmhost"]["compression"],
            bwa = "bwa-mem2" if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" else "bwa",
            minimum_seed_length = config["params"]["rmhost"]["bwa"]["minimum_seed_length"],
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            pe_bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}/pe"),
            se_bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}/se"),
            save_bam = "yes" if config["params"]["rmhost"]["save_bam"] else "no"
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["align"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR {params.report_dir}

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""
            STATSPE=""
            STATSSE=""

            if [ "$R1" != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz
                STATSPE={params.report_dir}/rmhost.pe.align_stats.txt

                mkdir -p $OUTPE
                rm -rf {params.pe_bam_dir}
                mkdir -p {params.pe_bam_dir}

                BAM={params.pe_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $R1 $R2 2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSPE) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 $FQ1 \
                        -2 $FQ2 -) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.pe_bam_dir}/temp \
                    -O BAM -o $BAM -
                else
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $R1 $R2 2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSPE) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 12 -F 256 \
                    -1 $FQ1 \
                    -2 $FQ2 -
                fi
            fi

            if [ "$RS" != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz
                STATSSE={params.report_dir}/rmhost.se.align_stats.txt

                mkdir -p $OUTSE
                rm -rf {params.se_bam_dir}
                mkdir -p {params.se_bam_dir}

                BAM={params.se_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $RS 2>> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSSE) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -cf -p {threads} \
                        > $FQS) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.se_bam_dir}/temp \
                    -O BAM -o $BAM -
                else
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $RS 2>> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSSE) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 4 -F 256 - | \
                    pigz -cf -p {threads} \
                    > $FQS
                fi
            fi

            echo "{{ \
            \\"PE_FORWARD\\": \\"$FQ1\\", \
            \\"PE_REVERSE\\": \\"$FQ2\\", \
            \\"SE\\": \\"$FQS\\", \
            \\"PE_ALIGN_STATS\\": \\"$STATSPE\\", \
            \\"SE_ALIGN_STATS\\": \\"$STATSSE\\" }}" | \
            jq . > {output}
            '''


    rule rmhost_bwa_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
                sample=SAMPLES_ID_LIST)

else:
    rule rmhost_bwa_all:
        input:


if config["params"]["rmhost"]["bowtie2"]["do"]:
    rule rmhost_bowtie2_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand(
                "{prefix}.{suffix}",
                prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_bowtie2_index/rmhost_bowtie2_index.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_bwotie2_index/rmhost_bowtie2_index.txt")
        params:
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"]
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["align"]
        shell:
            '''
            bowtie2-build --threads {threads} {input} {params.index_prefix} >{log} 2>&1
            '''


    rule rmhost_bowtie2:
        input:
            reads = lambda wildcards: rmhost_input(wildcards),
            index = expand(
                "{prefix}.{suffix}",
                prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_bowtie2/{sample}.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_bowtie2/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/se/{sample}"),
            report_dir = os.path.join(config["output"]["rmhost"], "report/flagstat/{sample}"),
            presets = config["params"]["rmhost"]["bowtie2"]["presets"],
            compression = config["params"]["rmhost"]["compression"],
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"],
            pe_bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}/pe"),
            se_bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}/se"),
            save_bam = "yes" if config["params"]["rmhost"]["save_bam"] else "no"
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["align"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR {params.report_dir}

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""
            STATSPE=""
            STATSSE=""

            if [ "$R1" != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz
                STATSPE={params.report_dir}/rmhost.pe.align_stats.txt

                mkdir -p $OUTPE
                rm -rf {params.pe_bam_dir}
                mkdir -p {params.pe_bam_dir}

                BAM={params.pe_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    bowtie2 \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -1 $R1 \
                    -2 $R2 \
                    {params.presets} \
                    2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - $STATSPE) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 $FQ1 \
                        -2 $FQ2 -) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.pe_bam_dir}/temp \
                    -O BAM -o $BAM -
                else
                    bowtie2 \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -1 $R1 \
                    -2 $R2 \
                    {params.presets} \
                    2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSPE) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 12 -F 256 \
                    -1 $FQ1 \
                    -2 $FQ2 -
                fi
            fi

            if [ "$RS" != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz
                STATSSE={params.report_dir}/rmhost.se.align_stats.txt

                mkdir -p $OUTSE
                rm -rf {params.se_bam_dir}
                mkdir -p {params.se_bam_dir}

                BAM={params.se_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    bowtie2 \
                    {params.presets} \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -U $RS \
                    2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSSE) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -cf -p {threads} \
                        > $FQS) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.se_bam_dir}/temp \
                    -O BAM -o $BAM -
                else
                    bowtie2 \
                    {params.presets} \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -U $RS \
                    2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSSE) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 4 -F 256 - | \
                    pigz -cf -p {threads} \
                    > $FQS
                fi
            fi

            echo "{{ \
            \\"PE_FORWARD\\": \\"$FQ1\\", \
            \\"PE_REVERSE\\": \\"$FQ2\\", \
            \\"SE\\": \\"$FQS\\", \
            \\"PE_ALIGN_STATS\\": \\"$STATSPE\\", \
            \\"SE_ALIGN_STATS\\": \\"$STATSSE\\" }}" | \
            jq . > {output}
            '''


    rule rmhost_bowtie2_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
                sample=SAMPLES_ID_LIST)

else:
    rule rmhost_bowtie2_all:
        input:


if config["params"]["rmhost"]["minimap2"]["do"]:
    rule rmhost_minimap2_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            config["params"]["rmhost"]["minimap2"]["index"]
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_minimap2_index/rmhost_minima2_index.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_minimap2_index/rmhost_minima2_index.txt")
        params:
            split_size = config["params"]["rmhost"]["minimap2"]["split_size"]
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["align"]
        shell:
            '''
            minimap2 -t {threads} -I {params.split_size} -d {output} {intput} >{log} 2>&1
            '''


    rule rmhost_minimap2:
        input:
            reads = lambda wildcards: rmhost_input(wildcards),
            index = config["params"]["rmhost"]["minimap2"]["index"]
        output:
            os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_minimap2/{sample}.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_minimap2/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/se/{sample}"),
            report_dir = os.path.join(config["output"]["rmhost"], "report/flagstat/{sample}"),
            preset = config["params"]["rmhost"]["minimap2"]["preset"],
            compression = config["params"]["rmhost"]["compression"],
            pe_bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}/pe"),
            se_bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}/se"),
            save_bam = "yes" if config["params"]["rmhost"]["save_bam"] else "no"
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["align"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR {params.report_dir}

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""
            STATSPE=""
            STATSSE=""

            if [ "$R1" != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz
                STATSPE={params.report_dir}/rmhost.pe.align_stats.txt

                mkdir -p $OUTPE
                rm -rf {params.pe_bam_dir}
                mkdir -p {params.pe_bam_dir}

                BAM={params.pe_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    minimap2 \
                    -t {threads} \
                    -ax {params.preset} \
                    {input.index} \
                    $R1 $R2 2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSPE) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 $FQ1 \
                        -2 $FQ2 -) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.pe_bam_dir}/temp \
                    -O BAM -o $BAM -
                else
                    minimap2 \
                    -t {threads} \
                    -ax {params.preset} \
                    {input.index} \
                    $R1 $R2 2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSPE) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 12 -F 256 \
                    -1 $FQ1 \
                    -2 $FQ2 -
                fi
            fi

            if [ "$RS" != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz
                STATSSE={params.report_dir}/rmhost.se.align_stats.txt

                mkdir -p $OUTSE
                rm -rf {params.se_bam_dir}
                mkdir -p {params.se_bam_dir}

                BAM={params.se_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    minimap2 \
                    -t {threads} \
                    {input.index} \
                    $RS 2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSSE) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -cf -p {threads} \
                        > $FQS) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.se_bam_dir}/temp \
                    -O BAM -o $BAM -
                else
                    minimap2 \
                    -t {threads} \
                    {input.index} \
                    $RS 2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - > $STATSSE) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 4 -F 256 - | \
                    pigz -cf -p {threads} \
                    > $FQS
                fi
            fi

            echo "{{ \
            \\"PE_FORWARD\\": \\"$FQ1\\", \
            \\"PE_REVERSE\\": \\"$FQ2\\", \
            \\"SE\\": \\"$FQS\\", \
            \\"PE_ALIGN_STATS\\": \\"$STATSPE\\", \
            \\"SE_ALIGN_STATS\\": \\"$STATSSE\\" }}" | \
            jq . > {output}
            '''


    rule rmhost_minimap2_all:
        input:
            expand(os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)

else:
    rule rmhost_minimap2_all:
        input:


if config["params"]["rmhost"]["kraken2"]["do"]:
    rule rmhost_kraken2:
        input:
            reads = lambda wildcards: rmhost_input(wildcards),
            database = config["params"]["rmhost"]["kraken2"]["database"]
        output:
            os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_kraken2/{sample}.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_kraken2/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/se/{sample}"),
            confidence = config["params"]["rmhost"]["kraken2"]["confidence"],
            min_base_quality = config["params"]["rmhost"]["kraken2"]["min_base_quality"],
            min_hit_groups = config["params"]["rmhost"]["kraken2"]["min_hit_groups"],
            host_taxid = config["params"]["rmhost"]["kraken2"]["host_taxid"]
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["kraken2"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""

            if [ "$R1" != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz
                TABLE={params.pe_prefix}.rmhost.pe.kraken2.table
                REPORT={params.pe_prefix}.rmhost.pe.kraken2.report.gz

                mkdir -p $OUTPE

                kraken2 \
                --threads {threads} \
                --db {input.database} \
                --use-names \
                --confidence {params.confidence} \
                --minimum-base-quality {params.min_base_quality} \
                --minimum-hit-groups {params.min_hit_groups} \
                --output $TABLE \
                --report $REPORT \
                --gzip-compressed \
                --paired \
                $R1 $R2 \
                >{log} 2>&1

                pigz -f -p {threads} ${{REPORT%.gz}}

                extract_kraken2_reads.py \
                -k $TABLE \
                --taxid {params.host_taxid} \
                --noappend \
                --exclude \
                --fastq-output \
                --gzip-output \
                -s $R1 \
                -s2 $R2 \
                -o $FQ1 \
                -o2 $FQ2 \
                >>{log} 2>&1
            fi

            if [ "$RS" != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz
                TABLE={params.se_prefix}.rmhost.se.kraken2.table
                REPORT={params.se_prefix}.rmhost.se.kraken2.report.gz

                mkdir -p $OUTSE

                kraken2 \
                --threads {threads} \
                --db {input.database} \
                --use-names \
                --confidence {params.confidence} \
                --minimum-base-quality {params.min_base_quality} \
                --minimum-hit-groups {params.min_hit_groups} \
                --output $TABLE \
                --report $REPORT \
                --gzip-compressed \
                $RS \
                >>{log} 2>&1

                pigz -f -p {threads} ${{REPORT%.gz}}

                extract_kraken2_reads.py \
                -k $TABLE \
                --taxid {params.host_taxid} \
                --noappend \
                --exclude \
                --fastq-output \
                --gzip-output \
                -s $RS \
                -o $FQS \
                >>{log} 2>&1
            fi

            echo "{{ \
            \\"PE_FORWARD\\": \\"$FQ1\\", \
            \\"PE_REVERSE\\": \\"$FQ2\\", \
            \\"SE\\": \\"$FQS\\ "}}" | \
            jq . > {output}
            '''


    rule rmhost_kraken2_all:
        input:
            expand(os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)

else:
    rule rmhost_kraken2_all:
        input:


if config["params"]["rmhost"]["kneaddata"]["do"]:
    rule rmhost_kneaddata:
        input:
            reads = lambda wildcards: rmhost_input(wildcards)
        output:
            os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json")
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_kneaddata/{sample}.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_kneaddata/{sample}.txt")
        params:
            pe_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/pe/{sample}"),
            se_prefix = os.path.join(config["output"]["rmhost"], "reads/{sample}/se/{sample}"),
            trf_options = "--run-trf" if config["params"]["rmhost"]["kneaddata"]["do_trf"] else "--bypass-trf",
            do_trimmomatic = "yes" if config["params"]["rmhost"]["kneaddata"]["do_trimmomatic"] else "no",
            trimmomatic_options = config["params"]["rmhost"]["kneaddata"]["trimmomatic_options"],
            sequencer_source = config["params"]["rmhost"]["kneaddata"]["sequencer_source"],
            do_bowtie2 = "yes" if config["params"]["rmhost"]["kneaddata"]["do_bowtie2"] else "no",
            bowtie2_options = config["params"]["rmhost"]["kneaddata"]["bowtie2_options"],
            decontaminate_pairs = config["params"]["rmhost"]["kneaddata"]["decontaminate_pairs"],
            bowtie2_database = config["params"]["rmhost"]["kneaddata"]["bowtie2_database"],
            do_bmtagger = "yes" if config["params"]["rmhost"]["kneaddata"]["do_bmtagger"] else "no",
            output_prefix = "{sample}.rmhost"
        priority:
            20
        threads:
            config["params"]["rmhost"]["threads"]
        conda:
            config["envs"]["kneaddata"]
        shell:
            '''
            OUTDIR=$(dirname {output})
            OUTPE=$(dirname {params.pe_prefix})
            OUTSE=$(dirname {params.se_prefix})

            rm -rf $OUTDIR $OUTPE $OUTSE
            mkdir -p $OUTDIR

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            FQ1=""
            FQ2=""
            FQS=""

            trimpath=$(which trimmomatic)
            trimpathreal=$(realpath $trimpath)
            trimdir=$(dirname $trimpathreal)

            if [ "$R1" != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz

                if [ "{params.do_bowtie2}" == "yes" ];
                then
                    if [ "{params.do_trimmomatic}" == "yes" ];
                    then
                        kneaddata \
                        -i $R1 -i $R2 \
                        {params.trf_options} \
                        --output $OUTPE \
                        --output-prefix {params.output_prefix} \
                        --reference-db {params.bowtie2_database} \
                        --trimmomatic $trimdir \
                        --trimmomatic-options '{params.trimmomatic_options}' \
                        --sequencer-source {params.sequencer_source} \
                        --bowtie2-options '{params.bowtie2_options} ' \
                        --decontaminate-pairs {params.decontaminate_pairs} \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                    else
                        kneaddata \
                        -i $R1 -i $R2 \
                        {params.trf_options} \
                        --bypass-trim \
                        --output $OUTPE \
                        --output-prefix {params.output_prefix} \
                        --reference-db {params.bowtie2_database} \
                        --bowtie2-options '{params.bowtie2_options} ' \
                        --decontaminate-pairs {params.decontaminate_pairs} \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                    fi
                else
                    if [ "{params.do_bmtagger}" == "yes" ];
                    then
                        if [ "{params.do_trimmomatic}" == "yes" ];
                        then
                            kneaddata \
                            -i $R1 -i $R2 \
                            {params.trf_options} \
                            --output $OUTPE \
                            --output-prefix {params.output_prefix} \
                            --trimmomatic $trimdir \
                            --trimmomatic-options '{params.trimmomatic_options}' \
                            --sequencer-source {params.sequencer_source} \
                            --run-bmtagger \
                            --remove-intermediate-output \
                            --threads {threads} \
                            --reorder \
                            --log {log}
                        else
                            kneaddata \
                            -i $R1 -i $R2 \
                            {params.trf_options} \
                            --bypass-trim \
                            --output $OUTPE \
                            --output-prefix {params.output_prefix} \
                            --run-bmtagger \
                            --remove-intermediate-output \
                            --threads {threads} \
                            --reorder \
                            --log {log}
                        fi
                    fi
                fi

                pigz -f -p {threads} $OUTPE/*
                mv $OUTPE/{params.output_prefix}_paired_1.fastq.gz $FQ1
                mv $OUTPE/{params.output_prefix}_paired_2.fastq.gz $FQ2
            fi

            if [ "$RS" != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz

                if [ "{params.do_bowtie2}" == "yes" ];
                then
                    if [ "{params.do_trimmomatic}" == "yes" ];
                    then
                        kneaddata \
                        -i $RS \
                        {params.trf_options} \
                        --output $OUTSE \
                        --output-prefix {params.output_prefix} \
                        --reference-db {params.bowtie2_database} \
                        --trimmomatic $trimdir \
                        --trimmomatic-options '{params.trimmomatic_options}' \
                        --sequencer-source {params.sequencer_source} \
                        --bowtie2-options '{params.bowtie2_options} ' \
                        --decontaminate-pairs {params.decontaminate_pairs} \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                    else
                        kneaddata \
                        -i $RS \
                        {params.trf_options} \
                        --bypass-trim \
                        --output $OUTSE \
                        --output-prefix {params.output_prefix} \
                        --reference-db {params.bowtie2_database} \
                        --bowtie2-options '{params.bowtie2_options} ' \
                        --decontaminate-pairs {params.decontaminate_pairs} \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                    fi
                else
                    if [ "{params.do_bmtagger}" == "yes" ];
                    then
                        if [ "{params.do_trimmomatic}" == "yes" ];
                        then
                            kneaddata \
                            -i $RS \
                            {params.trf_options} \
                            --output $OUTSE \
                            --output-prefix {params.output_prefix} \
                            --trimmomatic $trimdir \
                            --trimmomatic-options '{params.trimmomatic_options}' \
                            --sequencer-source {params.sequencer_source} \
                            --run-bmtagger \
                            --remove-intermediate-output \
                            --threads {threads} \
                            --reorder \
                            --log {log}
                        else
                            kneaddata \
                            -i $RS \
                            {params.trf_options} \
                            --bypass-trim \
                            --output $OUTSE \
                            --output-prefix {params.output_prefix} \
                            --run-bmtagger \
                            --remove-intermediate-output \
                            --threads {threads} \
                            --reorder \
                            --log {log}
                        fi
                    fi
                fi

                pigz -f -p {threads} $OUTSE/*
                mv $OUTSE/{params.output_prefix}.fastq.gz $FQS
            fi

            echo "{{ \
            \\"PE_FORWARD\\": \\"$FQ1\\", \
            \\"PE_REVERSE\\": \\"$FQ2\\", \
            \\"SE\\": \\"$FQS\\ "}}" | \
            jq . > {output}
            '''


    rule rmhost_kneaddata_all:
        input:
            expand(os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)

else:
    rule rmhost_kneaddata_all:
        input:


if RMHOST_DO \
and (not config["params"]["rmhost"]["kraken2"]["do"]) \
and (not config["params"]["rmhost"]["kneaddata"]["do"]):
    rule rmhost_alignment_report:
        input:
            expand(os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["rmhost"], "report/rmhost_align2host_stats.tsv")
        priority:
            20
        run:
            input_list = [str(i) for i in input]
            metapi.flagstats_summary(input_list, 2, output=str(output))


    rule rmhost_alignment_report_all:
        input:
            rules.rmhost_alignment_report.output

else:
    rule rmhost_alignment_report_all:
        input:


if RMHOST_DO and config["params"]["qcreport"]["do"]:
    rule rmhost_report:
        input:
            os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json")
        output:
            os.path.join(config["output"]["rmhost"], "report/stats/{sample}_rmhost_stats.tsv")
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_report/{sample}.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_report/{sample}.txt")
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


    rule rmhost_report_merge:
        input:
            expand(os.path.join(config["output"]["rmhost"], "report/stats/{sample}_rmhost_stats.tsv"),
            sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "rmhost_stats.tsv")
        log:
            os.path.join(config["output"]["raw"], "logs/rmhost_report_merge/rmhost_report_merge.log")
        benchmark:
            os.path.join(config["output"]["raw"], "benchmark/rmhost_report_merge/rmhost_report_merge.txt")
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


    rule rmhost_report_all:
        input:
            rules.rmhost_report_merge.output

else:
    rule rmhost_report_all:
        input:


rule rmhost_all:
    input:
        rules.rmhost_bwa_all.input,
        rules.rmhost_bowtie2_all.input,
        rules.rmhost_minimap2_all.input,
        rules.rmhost_kraken2_all.input,
        rules.rmhost_kneaddata_all.input,
        rules.rmhost_alignment_report_all.input,
        rules.rmhost_report_all.input


localrules:
    rmhost_bwa_all,
    rmhost_bowtie2_all,
    rmhost_minimap2_all,
    rmhost_kraken2_all,
    rmhost_kneaddata_all,
    rmhost_alignment_report_all,
    rmhost_report_all,
    rmhost_all