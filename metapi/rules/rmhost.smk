def rmhost_input(wildcards, have_single):
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming", have_single)
    else:
        return get_reads(wildcards, "raw", have_single)


if config["params"]["rmhost"]["bwa"]["do"]:
    rule rmhost_bwa_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                   suffix=["amb", "ann", "bwt", "pac", "sa"])
        log:
            os.path.join(config["output"]["rmhost"],
                         "logs/build_host_index_for_bwa.log")
        params:
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"]
        shell:
            '''
            bwa index {input} -p {params.index_prefix} 2> {log}
            '''


    rule rmhost_bwa:
        input:
            reads = lambda wildcards: rmhost_input(wildcards, False),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                           suffix=["amb", "ann", "bwt", "pac", "sa"])
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/flagstat/{sample}.align2host.flagstat"),
            reads = expand(os.path.join(config["output"]["rmhost"],
                                        "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "")
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.bwa.log")
        params:
            compression = config["params"]["rmhost"]["compression"],
            minimum_seed_length = config["params"]["rmhost"]["bwa"]["minimum_seed_length"],
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bam = os.path.join(config["output"]["rmhost"],
                               "bam/{sample}/{sample}.align2host.sorted.bam")
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["save_bam"]:
                    shell("mkdir -p %s" % os.path.dirname(params.bam))
                    shell(
                        '''
                        bwa mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} {input.reads[1]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 12 -F 256 \
                              -1 {output.reads[0]} \
                              -2 {output.reads[1]} -) | \
                        samtools sort \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} - \
                        2> {log}
                        ''')
                else:
                    shell(
                        '''
                        bwa mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} {input.reads[1]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 {output.reads} \
                        -2 {output.r2} - \
                        2> {log}
                        ''')
            else:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        bwa mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 4 -F 256 - \
                              > {output.reads[0]}) | \
                        samtools sort \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} - \
                        2> {log}
                        ''')
                else:
                    shell(
                        '''
                        bwa mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 4 -F 256 - \
                        > {output.reads[0]} \
                        2> {log}
                        ''')


    rule rmhost_bwa_all:
        input:
            expand([
                os.path.join(config["output"]["rmhost"],
                             "report/flagstat/{sample}.align2host.flagstat"),
                os.path.join(config["output"]["rmhost"],
                             "short_reads/{sample}/{sample}.rmhost{read}.fq.gz")],
                   read=[".1", ".2"] if IS_PE else "",
                   sample=SAMPLES.index.unique())


elif config["params"]["rmhost"]["bowtie2"]["do"]:
    rule rmhost_bowtie2_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                   suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        log:
            os.path.join(config["output"]["rmhost"], "logs/build_host_index_for_bowtie2.log")
        params:
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"]
        shell:
            '''
            bowtie2-build {input} {params.prefix} 2> {log}
            '''


    rule rmhost_bowtie2:
        input:
            reads = lambda wildcards: rmhost_input(wildcards, False),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/flagstat/{sample}.align2host.flagstat"),
            reads = expand(os.path.join(config["output"]["rmhost"],
                                        "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "")
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.bowtie2.log")
        params:
            compression = config["params"]["rmhost"]["compression"],
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"],
            bam = os.path.join(config["output"]["rmhost"],
                               "bam/{sample}/{sample}.align2host.sorted.bam")
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        bowtie2 \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -1 {input.reads[0]} \
                        -2 {input.reads[1]} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 12 -F 256 \
                              -1 {output.reads[0]} \
                              -2 {output.reads[1]} -) | \
                        samtools sort \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} -
                        ''')
                else:
                    shell(
                        '''
                        bowtie2 \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -1 {input.reads[0]} \
                        -2 {input.reads[1]} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 {output.reads[0]} \
                        -2 {output.reads[1]} -
                        ''')
            else:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        bowtie2 \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -U {input.reads[0]} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 4 -F 256 - \
                              > {output.reads[0]}) | \
                        samtools sort \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} -
                        ''')
                else:
                    shell(
                        '''
                        bowtie2 \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -U {input.reads[0]} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 4 -F 256 - \
                        > {output.reads[0]}
                        ''')


    rule rmhost_bowtie2_all:
        input:
            expand([
                os.path.join(config["output"]["rmhost"],
                             "report/flagstat/{sample}.align2host.flagstat"),
                os.path.join(config["output"]["rmhost"],
                             "short_reads/{sample}/{sample}.rmhost{read}.fq.gz")],
                   read=[".1", ".2"] if IS_PE else "",
                   sample=SAMPLES.index.unique())


if RMHOST_DO:
    rule rmhost_report:
        input:
            expand(os.path.join(config["output"]["rmhost"],
                                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                   read=[".1", ".2"] if IS_PE else "")
        output:
            os.path.join(config["output"]["rmhost"],
                         "report/stats/{sample}_rmhost_stats.tsv")
        params:
            fq_encoding = config["params"]["fq_encoding"],
            sample_id = "{sample}"
        threads:
            config["params"]["qc_report"]["seqkit"]["threads"]
        run:
            if IS_PE:
                shell(
                    '''
                    seqkit stats \
                    --all --basename --tabular \
                    --fq-encoding %s \
                    --out-file %s \
                    --threads %d %s
                    ''' % (params.fq_encoding, output, threads, " ".join(input)))
                metapi.change(output[0], params.sample_id, "rmhost", "pe", ["fq1", "fq2"])
            else:
                shell(
                    '''
                    seqkit stats \
                    --all --basename --tabular \
                    --fq-encoding %s \
                    --out-file %s \
                    --threads %d %s
                    ''' % (params.fq_encoding, output, threads, input))
                metapi.change(output[0], params.sample_id, "rmhost", "se", ["fq1"])


    rule rmhost_report_merge:
        input:
            expand(os.path.join(config["output"]["rmhost"],
                                "report/stats/{sample}_rmhost_stats.tsv"),
                   sample=SAMPLES.index.unique())
        output:
            os.path.join(config["output"]["rmhost"],
                         "report/rmhost_stats.tsv")
        threads:
            config["params"]["qc_report"]["seqkit"]["threads"]
        run:
            metapi.merge(input, metapi.parse, threads, save=True, output=output[0])


    rule rmhost_report_all:
        input:
            os.path.join(config["output"]["rmhost"], "report/rmhost_stats.tsv")


rule qc_report:
    input:
        #raw_stats = os.path.join(config["results"]["report"]["base_dir"], "raw.stats.tsv"),
        trim_stats = os.path.join(config["output"]["trimming"],
                                  "report/trimming_stats.tsv"),
        rmhost_stats = os.path.join(config["output"]["rmhost"],
                                    "report/rmhost_stats.tsv")
    output:
        os.path.join(config["output"]["rmhost"],
                     "report/qc_stats.tsv")
    threads:
        config["params"]["qc_report"]["seqkit"]["threads"]
    run:
        df = metapi.merge([input.trim_stats, input.rmhost_stats],
                          metapi.parse, threads)
        metapi.compute_host_rate(df, save=True, output=output[0])


rule qc_report_all:
    input:
        os.path.join(config["output"]["rmhost"],
                     "report/qc_stats.tsv")


if config["params"]["rmhost"]["bwa"]["do"]:
    rule rmhost_all:
        input:
            rules.rmhost_bwa_all.input,
            rules.rmhost_report_all.input,
            rules.qc_report_all.input

elif config["params"]["rmhost"]["bowtie2"]["do"]:
    rule rmhost_all:
        input:
            rules.rmhost_bowtie2_all.input,
            rules.rmhost_report_all.input,
            rules.qc_report_all.input

else:
    rule rmhost_all:
        input:
