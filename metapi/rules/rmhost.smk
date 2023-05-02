def rmhost_input(wildcards):
    if TRIMMING_DO:
        return os.path.join(config["output"]["trimming"], f"reads/{wildcards.sample}/{wildcards.sample}.json")
    else:
        return os.path.join(config["output"]["raw"], f"reads/{wildcards.sample}/{wildcards.sample}.json")


BWA_INDEX_SUFFIX = ["0123", "amb", "ann", "bwt.2bit.64", "pac"] if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" \
else ["amb", "ann", "bwt", "pac", "sa"]

if config["params"]["rmhost"]["bwa"]["do"]:
    rule rmhost_bwa_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}", prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
            suffix=BWA_INDEX_SUFFIX)
        log:
            os.path.join(config["output"]["rmhost"], "logs/build_host_index_for_bwa.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/bwa.index.benchmark.txt")
        params:
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bwa = "bwa-mem2" if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" else "bwa"
        priority:
            20
        conda:
            config["envs"]["align"]
        shell:
            '''
            {params.bwa} index -p {params.index_prefix} {input} 2> {log}
            '''


    rule rmhost_bwa:
        input:
            reads = lambda wildcards: rmhost_input(wildcards),
            index = expand("{prefix}.{suffix}", prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
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
            mkdir -p $OUTDIR $OUTPE $OUTSE

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            if [ $R1 != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz
                STATS={params.pe_prefix}.pe.align_stats.txt

                rm -rf {params.pe_bam_dir}
                mkdir -p {params.pe_bam_dir}

                BAM={params.pe_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $R1 $R2 | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > $STATS) | \
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
                    -O BAM -o $BAM - \
                    2> {log}
                else
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $R1 $R2 | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > $STATS) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 12 -F 256 \
                    -1 $FQ1 \
                    -2 $FQ2 - \
                    2> {log}
                fi
            fi

            if [ $RS != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz
                STATS={params.se_prefix}.se.align_stats.txt

                rm -rf {params.se_bam_dir}
                mkdir -p {params.se_bam_dir}

                BAM={params.se_bam_dir}/sorted.bam

                if [ "{params.save_bam}" == "yes" ];
                then
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $RS | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > $STATS) | \
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
                    -O BAM -o $BAM - \
                    2>> {log}
                else
                    {params.bwa} mem \
                    -k {params.minimum_seed_length} \
                    -t {threads} \
                    {params.index_prefix} \
                    $RS | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > $STATS) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 4 -F 256 - | \
                    pigz -cf -p {threads} \
                    > $FQS \
                    2>> {log}
                fi
            fi

            echo "{{\\"PE_FORWARD\\": \\"$FQ1\\", \\"PE_REVERSE\\": \\"$FQ2\\", \\"SE\\": \\"$FQS\\"}}" | \
            jq . > {output}
            '''


    rule rmhost_bwa_all:
        input:
            expand(os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
            sample=SAMPLES_ID_LIST)

else:
    rule rmhost_bwa_all:
        input:


if config["params"]["rmhost"]["bowtie2"]["do"]:
    rule rmhost_bowtie2_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}", prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        log:
            os.path.join(config["output"]["rmhost"], "logs/rmhost_bowtie2_index/rmhost_bowtie2_index.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/rmhost_bowtie2_index.txt")
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
            index = expand("{prefix}.{suffix}", prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
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
            mkdir -p $OUTDIR $OUTPE $OUTSE

            R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
            R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
            RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

            if [ $R1 != "" ];
            then
                FQ1={params.pe_prefix}.rmhost.pe.1.fq.gz
                FQ2={params.pe_prefix}.rmhost.pe.2.fq.gz
                STATS={params.pe_prefix}.pe.align_stats.txt

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
                        -@4 - | \
                        pigz -cf > $STATS) | \
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
                        -@4 - | \
                        pigz -cf > $STATS) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 12 -F 256 \
                    -1 $FQ1 \
                    -2 $FQ2 -
                fi
            fi

            if [ $RS != "" ];
            then
                FQS={params.se_prefix}.rmhost.se.fq.gz
                STATS={params.se_prefix}.se.align_stats.txt

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
                        -@4 - | \
                        pigz -cf > $STATS) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -cf -p {threads} \
                        > $FQS) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.se_bam_dir} \
                    -O BAM -o $BAM -
                else
                    bowtie2 \
                    {params.presets} \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -U $RS \
                    2> {log} | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > $STATS) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 4 -F 256 - | \
                    pigz -cf -p {threads} \
                    > $FQS
                fi
            fi

            echo "{{\\"PE_FORWARD\\": \\"$FQ1\\", \\"PE_REVERSE\\": \\"$FQ2\\", \\"SE\\": \\"$FQS\\"}}" | \
            jq . > {output}
            '''


    rule rmhost_bowtie2_all:
        input:
            expand(os.path.join(config["output"]["rmhost"], "reads/{sample}/{sample}.json"),
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
        params:
            split_size = config["params"]["rmhost"]["minimap2"]["split_size"]
        conda:
            config["envs"]["align"]
        priority:
            20
        shell:
            '''
            minimap2 -I {params.split_size} -d {output} {intput}
            '''


    if IS_PE:
        rule rmhost_minimap2:
            input:
                lambda wildcards: trimming_stats_input(wildcards),
                reads = lambda wildcards: rmhost_input(wildcards),
                index = config["params"]["rmhost"]["minimap2"]["index"]
            output:
                flagstat = os.path.join(config["output"]["rmhost"],
                                        "report/flagstat/{sample}.align2host.flagstat.gz"),
                reads = expand(os.path.join(
                    config["output"]["rmhost"],
                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                            read=[".1", ".2"]) \
                            if config["params"]["rmhost"]["save_reads"] else \
                                temp(expand(os.path.join(
                                    config["output"]["rmhost"],
                                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                            read=[".1", ".2"]))
            conda:
                config["envs"]["align"]
            log:
                os.path.join(config["output"]["rmhost"], "logs/{sample}.minimap2.log")
            benchmark:
                os.path.join(config["output"]["rmhost"],
                            "benchmark/minimap2/{sample}.minimap2.txt")
            priority:
                20
            params:
                preset = config["params"]["rmhost"]["minimap2"]["preset"],
                compression = config["params"]["rmhost"]["compression"],
                bam = os.path.join(config["output"]["rmhost"], "bam/{sample}/{sample}.align2host.sorted.bam"),
                bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}"),
                save_bam = "yes" if config["params"]["rmhost"]["save_bam"] else "no"
            threads:
                config["params"]["rmhost"]["threads"]
            shell:
                '''
                if [ "{params.save_bam}" == "yes" ];
                then
                    mkdir -p {params.bam_dir}

                    minimap2 \
                    -t {threads} \
                    -ax {params.preset} \
                    {input.index} \
                    {input.reads[0]} {input.reads[1]} | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > {output.flagstat}) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 {output.reads[0]} \
                        -2 {output.reads[1]} -) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.bam} \
                    -O BAM -o {params.bam} - \
                    2> {log}
                else
                    minimap2 \
                    -t {threads} \
                    -ax {params.preset} \
                    {input.index} \
                    {input.reads[0]} {input.reads[1]} | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > {output.flagstat}) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 12 -F 256 \
                    -1 {output.reads[0]} \
                    -2 {output.reads[1]} - \
                    2> {log}
                fi
                '''

    else:
        rule rmhost_minimap2:
            input:
                lambda wildcards: trimming_stats_input(wildcards),
                reads = lambda wildcards: rmhost_input(wildcards),
                index = config["params"]["rmhost"]["minimap2"]["index"]
            output:
                flagstat = os.path.join(config["output"]["rmhost"],
                                        "report/flagstat/{sample}.align2host.flagstat.gz"),
                reads = expand(os.path.join(
                    config["output"]["rmhost"],
                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                            read=[""]) \
                            if config["params"]["rmhost"]["save_reads"] else \
                                temp(expand(os.path.join(
                                    config["output"]["rmhost"],
                                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                            read=[""]))
            conda:
                config["envs"]["align"]
            log:
                os.path.join(config["output"]["rmhost"], "logs/{sample}.minimap2.log")
            benchmark:
                os.path.join(config["output"]["rmhost"],
                            "benchmark/minimap2/{sample}.minimap2.txt")
            priority:
                20
            params:
                preset = config["params"]["rmhost"]["minimap2"]["preset"],
                compression = config["params"]["rmhost"]["compression"],
                bam = os.path.join(config["output"]["rmhost"], "bam/{sample}/{sample}.align2host.sorted.bam"),
                bam_dir = os.path.join(config["output"]["rmhost"], "bam/{sample}"),
                save_bam = "yes" if config["params"]["rmhost"]["save_bam"] else "no"
            threads:
                config["params"]["rmhost"]["threads"]
            shell:
                '''
                if [ "{params.save_bam}" == "yes" ];
                then
                    mkdir -p {params.bam_dir}

                    minimap2 \
                    -t {threads} \
                    {input.index} \
                    {input.reads[0]} | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > {output.flagstat}) | \
                    tee >(samtools fastq \
                        -@4 \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -cf -p {threads} \
                        > {output.reads[0]}) | \
                    samtools sort \
                    -m 3G \
                    -@4 \
                    -T {params.bam} \
                    -O BAM -o {params.bam} - \
                    2> {log}
                else
                    minimap2 \
                    -t {threads} \
                    {input.index} \
                    {input.reads[0]} | \
                    tee >(samtools flagstat \
                        -@4 - | \
                        pigz -cf > {output.flagstat}) | \
                    samtools fastq \
                    -@4 \
                    -c {params.compression} \
                    -N -f 4 -F 256 - | \
                    pigz -cf -p {threads} \
                    > {output.reads[0]} \
                    2> {log}
                fi
                '''


    rule rmhost_minimap2_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_minimap2_all:
        input:


if config["params"]["rmhost"]["kraken2"]["do"]:
    if IS_PE:
        rule rmhost_kraken2:
            input:
                lambda wildcards: trimming_stats_input(wildcards),
                reads = lambda wildcards: rmhost_input(wildcards),
                database = config["params"]["rmhost"]["kraken2"]["database"]
            output:
                table = temp(os.path.join(config["output"]["rmhost"],
                                        "short_reads/{sample}/{sample}.kraken2.table")),
                report = os.path.join(config["output"]["rmhost"],
                                    "short_reads/{sample}/{sample}.kraken2.report.gz"),
                reads = expand(os.path.join(
                    config["output"]["rmhost"],
                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                            read=[".1", ".2"]) \
                            if config["params"]["rmhost"]["save_reads"] else \
                                temp(expand(os.path.join(
                                    config["output"]["rmhost"],
                                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                            read=[".1", ".2"]))
            conda:
                config["envs"]["kraken2"]
            log:
                os.path.join(config["output"]["rmhost"], "logs/{sample}.kraken2.log")
            benchmark:
                os.path.join(config["output"]["rmhost"],
                            "benchmark/kraken2/{sample}.kraken2.txt")
            params:
                report = os.path.join(config["output"]["rmhost"],
                                    "short_reads/{sample}/{sample}.kraken2.report"),
                confidence = config["params"]["rmhost"]["kraken2"]["confidence"],
                min_base_quality = config["params"]["rmhost"]["kraken2"]["min_base_quality"],
                min_hit_groups = config["params"]["rmhost"]["kraken2"]["min_hit_groups"],
                host_taxid = config["params"]["rmhost"]["kraken2"]["host_taxid"]
            priority:
                20
            threads:
                config["params"]["rmhost"]["threads"]
            shell:
                '''
                kraken2 \
                --threads {threads} \
                --db {input.database} \
                --use-names \
                --confidence {params.confidence} \
                --minimum-base-quality {params.min_base_quality} \
                --minimum-hit-groups {params.min_hit_groups} \
                --output {output.table} \
                --report {params.report} \
                --gzip-compressed \
                --paired \
                {input.reads} \
                >{log} 2>&1

                pigz -f -p {threads} {params.report}

                extract_kraken2_reads.py \
                -k {output.table} \
                --taxid {params.host_taxid} \
                --noappend \
                --exclude \
                --fastq-output \
                --gzip-output \
                -s {input.reads[0]} \
                -s2 {input.reads[1]} \
                -o {output.reads[0]} \
                -o2 {output.reads[1]} \
                >>{log} 2>&1
                '''

    else:
        rule rmhost_kraken2:
            input:
                lambda wildcards: trimming_stats_input(wildcards),
                reads = lambda wildcards: rmhost_input(wildcards),
                database = config["params"]["rmhost"]["kraken2"]["database"]
            output:
                table = temp(os.path.join(config["output"]["rmhost"],
                                        "short_reads/{sample}/{sample}.kraken2.table")),
                report = os.path.join(config["output"]["rmhost"],
                                    "short_reads/{sample}/{sample}.kraken2.report.gz"),
                reads = expand(os.path.join(
                    config["output"]["rmhost"],
                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                            read=[""]) \
                            if config["params"]["rmhost"]["save_reads"] else \
                                temp(expand(os.path.join(
                                    config["output"]["rmhost"],
                                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                            read=[""]))
            conda:
                config["envs"]["kraken2"]
            log:
                os.path.join(config["output"]["rmhost"], "logs/{sample}.kraken2.log")
            benchmark:
                os.path.join(config["output"]["rmhost"],
                            "benchmark/kraken2/{sample}.kraken2.txt")
            params:
                report = os.path.join(config["output"]["rmhost"],
                                    "short_reads/{sample}/{sample}.kraken2.report"),
                confidence = config["params"]["rmhost"]["kraken2"]["confidence"],
                min_base_quality = config["params"]["rmhost"]["kraken2"]["min_base_quality"],
                min_hit_groups = config["params"]["rmhost"]["kraken2"]["min_hit_groups"],
                host_taxid = config["params"]["rmhost"]["kraken2"]["host_taxid"]
            priority:
                20
            threads:
                config["params"]["rmhost"]["threads"]
            shell:
                '''
                kraken2 \
                --threads {threads} \
                --db {input.database} \
                --use-names \
                --confidence {params.confidence} \
                --minimum-base-quality {params.min_base_quality} \
                --minimum-hit-groups {params.min_hit_groups} \
                --output {output.table} \
                --report {params.report} \
                --gzip-compressed \
                {input.reads[0]} \
                >{log} 2>&1

                pigz -f -p {threads} {params.report}

                extract_kraken2_reads.py \
                -k {output.table} \
                --taxid {params.host_taxid} \
                --noappend \
                --exclude \
                --fastq-output \
                --gzip-output \
                -s {input.reads[0]} \
                -o {output.reads[0]} \
                >>{log} 2>&1
                '''


    rule rmhost_kraken2_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_kraken2_all:
        input:


if config["params"]["rmhost"]["kneaddata"]["do"]:
    if IS_PE:
        rule rmhost_kneaddata:
            input:
                lambda wildcards: trimming_stats_input(wildcards),
                reads = lambda wildcards: rmhost_input(wildcards)
            output:
                reads = expand(os.path.join(
                    config["output"]["rmhost"],
                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                            read=[".1", ".2"]) \
                            if config["params"]["rmhost"]["save_reads"] else \
                                temp(expand(os.path.join(
                                    config["output"]["rmhost"],
                                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                            read=[".1", ".2"]))
            conda:
                config["envs"]["kneaddata"]
            log:
                os.path.join(config["output"]["rmhost"], "logs/{sample}.kneaddata.log")
            benchmark:
                os.path.join(config["output"]["rmhost"],
                            "benchmark/kneaddata/{sample}.kneaddata.txt")
            params:
                trf_options = "--run-trf" if config["params"]["rmhost"]["kneaddata"]["do_trf"] else "--bypass-trf",
                do_trimmomatic = "yes" if config["params"]["rmhost"]["kneaddata"]["do_trimmomatic"] else "no",
                trimmomatic_options = config["params"]["rmhost"]["kneaddata"]["trimmomatic_options"],
                sequencer_source = config["params"]["rmhost"]["kneaddata"]["sequencer_source"],
                do_bowtie2 = "yes" if config["params"]["rmhost"]["kneaddata"]["do_bowtie2"] else "no",
                bowtie2_options = config["params"]["rmhost"]["kneaddata"]["bowtie2_options"],
                decontaminate_pairs = config["params"]["rmhost"]["kneaddata"]["decontaminate_pairs"],
                bowtie2_database = config["params"]["rmhost"]["kneaddata"]["bowtie2_database"],
                do_bmtagger = "yes" if config["params"]["rmhost"]["kneaddata"]["do_bmtagger"] else "no",
                output_dir = os.path.join(config["output"]["rmhost"], "short_reads/{sample}"),
                output_prefix = "{sample}.rmhost",
            priority:
                20
            threads:
                config["params"]["rmhost"]["threads"]
            shell:
                '''
                rm -rf {params.output_dir}

                input_reads="-i {input.reads[0]} -i {input.reads[1]}"

                if [ "{params.do_bowtie2}" == "yes" ];
                then
                    if [ "{params.do_trimmomatic}" == "yes" ];
                    then
                        trimpath=$(which trimmomatic)
                        trimpathreal=$(realpath $trimpath)
                        trimdir=$(dirname $trimpathreal)

                        kneaddata $input_reads \
                        {params.trf_options} \
                        --output {params.output_dir} \
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
                        kneaddata $input_reads \
                        {params.trf_options} \
                        --bypass-trim \
                        --output {params.output_dir} \
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
                            trimpath=$(which trimmomatic)
                            trimpathreal=$(realpath $trimpath)
                            trimdir=$(dirname $trimpathreal)

                            kneaddata $input_reads \
                            {params.trf_options} \
                            --output {params.output_dir} \
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
                            kneaddata $input_reads \
                            {params.trf_options} \
                            --bypass-trim \
                            --output {params.output_dir} \
                            --output-prefix {params.output_prefix} \
                            --run-bmtagger \
                            --remove-intermediate-output \
                            --threads {threads} \
                            --reorder \
                            --log {log}
                        fi
                    fi
                fi

                pigz -f -p {threads} {params.output_dir}/*


                mv {params.output_dir}/{params.output_prefix}_paired_1.fastq.gz {output.reads[0]}
                mv {params.output_dir}/{params.output_prefix}_paired_2.fastq.gz {output.reads[1]}
                '''

    else:
        rule rmhost_kneaddata:
            input:
                lambda wildcards: trimming_stats_input(wildcards),
                reads = lambda wildcards: rmhost_input(wildcards)
            output:
                reads = expand(os.path.join(
                    config["output"]["rmhost"],
                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                            read=[""]) \
                            if config["params"]["rmhost"]["save_reads"] else \
                                temp(expand(os.path.join(
                                    config["output"]["rmhost"],
                                    "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                            read=[""]))
            conda:
                config["envs"]["kneaddata"]
            log:
                os.path.join(config["output"]["rmhost"], "logs/{sample}.kneaddata.log")
            benchmark:
                os.path.join(config["output"]["rmhost"],
                            "benchmark/kneaddata/{sample}.kneaddata.txt")
            params:
                trf_options = "--run-trf" if config["params"]["rmhost"]["kneaddata"]["do_trf"] else "--bypass-trf",
                do_trimmomatic = "yes" if config["params"]["rmhost"]["kneaddata"]["do_trimmomatic"] else "no",
                trimmomatic_options = config["params"]["rmhost"]["kneaddata"]["trimmomatic_options"],
                sequencer_source = config["params"]["rmhost"]["kneaddata"]["sequencer_source"],
                do_bowtie2 = "yes" if config["params"]["rmhost"]["kneaddata"]["do_bowtie2"] else "no",
                bowtie2_options = config["params"]["rmhost"]["kneaddata"]["bowtie2_options"],
                decontaminate_pairs = config["params"]["rmhost"]["kneaddata"]["decontaminate_pairs"],
                bowtie2_database = config["params"]["rmhost"]["kneaddata"]["bowtie2_database"],
                do_bmtagger = "yes" if config["params"]["rmhost"]["kneaddata"]["do_bmtagger"] else "no",
                output_dir = os.path.join(config["output"]["rmhost"], "short_reads/{sample}"),
                output_prefix = "{sample}.rmhost",
            priority:
                20
            threads:
                config["params"]["rmhost"]["threads"]
            shell:
                '''
                rm -rf {params.output_dir}

                input_reads="-i {input.reads}"

                if [ "{params.do_bowtie2}" == "yes" ];
                then
                    if [ "{params.do_trimmomatic}" == "yes" ];
                    then
                        trimpath=$(which trimmomatic)
                        trimpathreal=$(realpath $trimpath)
                        trimdir=$(dirname $trimpathreal)

                        kneaddata $input_reads \
                        {params.trf_options} \
                        --output {params.output_dir} \
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
                        kneaddata $input_reads \
                        {params.trf_options} \
                        --bypass-trim \
                        --output {params.output_dir} \
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
                            trimpath=$(which trimmomatic)
                            trimpathreal=$(realpath $trimpath)
                            trimdir=$(dirname $trimpathreal)

                            kneaddata $input_reads \
                            {params.trf_options} \
                            --output {params.output_dir} \
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
                            kneaddata $input_reads \
                            {params.trf_options} \
                            --bypass-trim \
                            --output {params.output_dir} \
                            --output-prefix {params.output_prefix} \
                            --run-bmtagger \
                            --remove-intermediate-output \
                            --threads {threads} \
                            --reorder \
                            --log {log}
                        fi
                    fi
                fi

                pigz -f -p {threads} {params.output_dir}/*


                mv {params.output_dir}/{params.output_prefix}.fastq.gz {output.reads}
                '''


    rule rmhost_kneaddata_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_kneaddata_all:
        input:


if RMHOST_DO \
and (not config["params"]["rmhost"]["kraken2"]["do"]) \
and (not config["params"]["rmhost"]["kneaddata"]["do"]):
    rule rmhost_alignment_report:
        input:
            expand(
                os.path.join(config["output"]["rmhost"],
                             "report/flagstat/{sample}.align2host.flagstat.gz"),
                sample=SAMPLES_ID_LIST)
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/rmhost_align2host_stats.tsv.gz")
        priority:
            20
        run:
            input_list = [str(i) for i in input]
            metapi.flagstats_summary(input_list, 2, output=output.flagstat)

else:
    rule rmhost_alignment_report:
        input:


if RMHOST_DO and config["params"]["qcreport"]["do"]:
    rule rmhost_report:
        input:
            lambda wildcards: get_reads(wildcards, "rmhost")
        output:
            temp(os.path.join(config["output"]["rmhost"],
                              "report/stats/{sample}_rmhost_stats.tsv.raw"))
        conda:
            config["envs"]["report"]
        log:
            os.path.join(config["output"]["rmhost"],
                         "logs/{sample}.seqkit.log")
        priority:
            25
        params:
            fq_encoding = config["params"]["fq_encoding"]
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        shell:
            '''
            seqkit stats \
            --all \
            --basename \
            --tabular \
            --fq-encoding {params.fq_encoding} \
            --out-file {output} \
            --threads {threads} \
            {input} 2> {log}
            '''


    rule rmhost_report_refine:
        input:
            os.path.join(config["output"]["rmhost"],
                         "report/stats/{sample}_rmhost_stats.tsv.raw")
        output:
            os.path.join(config["output"]["rmhost"],
                         "report/stats/{sample}_rmhost_stats.tsv.gz")
        params:
            sample_id = "{sample}"
        threads:
            1
        run:
            if IS_PE:
                metapi.change(str(input), str(output), params.sample_id, "rmhost",
                              "pe", ["fq1", "fq2"])
            else:
                metapi.change(str(input), str(output), params.sample_id, "rmhost",
                              "se", ["fq1"])


    rule rmhost_report_merge:
        input:
            expand(
                os.path.join(config["output"]["rmhost"],
                             "report/stats/{sample}_rmhost_stats.tsv.gz"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "rmhost_stats.tsv.gz")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            metapi.merge(input, metapi.parse, threads, output=output[0])


    rule rmhost_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "rmhost_stats.tsv.gz")

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
        rules.rmhost_alignment_report.input,
        rules.rmhost_report_all.input


localrules:
    rmhost_bwa_all,
    rmhost_bowtie2_all,
    rmhost_kraken2_all,
    rmhost_minimap2_all,
    rmhost_kneaddata_all,
    rmhost_report_all,
    rmhost_all