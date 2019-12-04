def trimmed_reads(wildcards, have_single):
    if have_single:
        return temp(expand(os.path.join(config["results"]["trimming"], "{sample}.trimmed{read}.fq.gz"),
                           sample=wildcards.sample,
                           read=[".1", ".2", ".single"] if IS_PE else ""))
    else:
        return temp(expand(os.path.join(config["results"]["trimming"], "{sample}.trimmed{read}.fq.gz"),
                           sample=wildcards.sample,
                           read=[".1", ".2"] if IS_PE else ""))


if config["params"]["rmhost"]["bwa"]["do"]:
    rule build_host_index_for_bwa:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                   suffix=["amb", "ann", "bwt", "pac", "sa"])
        log:
            os.path.join(config["logs"]["rmhost"], "build_host_index_for_bwa.log")
        params:
            prefix = config["params"]["rmhost"]["bwa"]["index_prefix"]
        shell:
            '''
            bwa index {input} -p {params.prefix} >{log} 2>&1
            '''

    rule rmhost_bwa:
        input:
            reads = lambda wildcards: trimmed_reads(wildcards, False),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                           suffix=["amb", "ann", "bwt", "pac", "sa"])
        output:
            flagstat = os.path.join(config["results"]["rmhost"], "{sample}.rmhost.flagstat.txt"),
            reads = expand(os.path.join(config["results"]["rmhost"], "{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "")
        log:
            os.path.join(config["logs"]["rmhost"], "{sample}.bwa.rmhost.log")
        params:
            minimum_seed_length = config["params"]["rmhost"]["bwa"]["minimum_seed_length"],
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bam = os.path.join(config["results"]["rmhost"], "{sample}.bwa.host.sorted.bam")
        threads:
            config["params"]["rmhost"]["bwa"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["bwa"]["save_bam"]:
                    shell('''bwa mem -k {params.minimum_seed_length} -t {threads} {params.index_prefix} {input.reads[0]} {input.reads[1]} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          tee >(samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.reads[0]} -2 {output.reads[1]} -) | \
                          samtools sort -@{threads} -O BAM -o {params.bam} - 2>{log}''')
                else:
                    shell('''bwa mem -k {params.minimum_seed_length} -t {threads} {params.prefix} {input.reads[0]} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.reads} -2 {output.r2} - 2>{log}''')
            else:
                if config["params"]["rmhost"]["bwa"]["save_bam"]:
                    shell('''bwa mem -k {params.minimum_seed_length} -t {threads} {params.index_prefix} {input.reads[0] | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat} | \
                          tee >(samtools fastq -@{threads} -N -f 8 -F 256 - > {output.reads[0]}) | \
                          samtools sort -@{threads} -O BAM -o {params.bam} - 2>{log}''')
                else:
                    shell('''bwa mem -k {params.minimum_seed_length} -t {threads} {params.prefix} {input.reads[0]} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          samtools fastq -@{threads} -N -f 8 -F 256 - > {output.reads[0]} 2>{log}''')

if config["params"]["rmhost"]["bowtie2"]["do"]:
    rule build_host_index_for_bowtie2:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                   suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        log:
            os.path.join(config["logs"]["rmhost"], "build_host_index_for_bowtie2.log")
        params:
            prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"]
        shell:
            '''
            bowtie2-build {input} {params.prefix} >{log} 2>&1
            '''

    rule rmhost_bowtie2:
        input:
            reads = lambda wildcards: trimmed_reads(wildcards, False),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            flagstat = os.path.join(config["results"]["rmhost"], "{sample}.rmhost.flagstat.txt"),
            reads = expand(os.path.join(config["results"]["rmhost"], "{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "")
        log:
            os.path.join(config["logs"]["rmhost"], "{sample}.bowtie2.rmhost.log")
        params:
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"],
            additional_params = config["params"]["rmhost"]["bowtie2"]["additional_params"],
            bam = os.path.join(config["results"]["rmhost"], "{sample}.bowtie2.host.sorted.bam")
        threads:
            config["params"]["rmhost"]["bowtie2"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["bowtie2"]["save_bam"]:
                    shell('''bowtie2 --threads {threads} -x {params.index_prefix} \
                          -1 {input.reads[0]} -2 {input.reads[1]} {params.additional_params} 2> {log} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          tee >(samtools sort -@{threads} -O BAM -o {params.bam}) | \
                          samtools view -@{threads} -SF4 - | awk -F'[/\t]' '{{print $1}}' | sort | uniq | \
                          tee >(awk '{{print $0 "/1"}}' - | seqtk subseq -r {input.r1} - | pigz -p {threads} -c > {output.reads[0]}) | \
                          awk '{{print $0 "/2"}}' - | seqtk subseq -r {input.r2} - | pigz -p {threads} -c > {output.reads[1]}''')
                else:
                    shell('''bowtie2 --threads {threads} -x {params.index_prefix} \
                          -1 {input.reads[0]} -2 {input.reads[1]} {params.additional_params} 2> {log} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          samtools view -@{threads} -SF4 - | awk -F'[/\t]' '{{print $1}}' | sort | uniq | \
                          tee >(awk '{{print $0 "/1"}}' - | seqtk subseq -r {input.reads[0]} - | pigz -p {threads} -c > {output.reads[0]}) | \
                          awk '{{print $0 "/2"}}' - | seqtk subseq -r {input.reads[1]} - | pigz -p {threads} -c > {output.reads[1]}''')
            else:
                if config["params"]["rmhost"]["bowtie2"]["save_bam"]:
                    shell('''bowtie2 --threads {threads} -x {params.index_prefix} \
                          -U {input.reads[0]} {params.additional_params} 2> {log} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          tee >(samtools sort -@{threads} -O BAM -o {params.bam}) | \
                          samtools view -@{threads} -SF4 - | awk -F'[/\t]' '{{print $1}}' | sort | uniq | \
                          seqtk subseq -r {input.reads[0]} - | pigz -p {threads} -c > {output.reads[0]}''')
                else:
                    shell('''bowtie2 --threads {threads} -x {params.index_prefix} \
                          -U {input.reads[0]} {params.additional_params} 2> {log} | \
                          tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                          samtools view -@{threads} -SF4 - | awk -F'[/\t]' '{{print $1}}' | sort | uniq | \
                          seqtk subseq -r {input.reads[0]} - | pigz -p {threads} -c > {output.reads[0]}''')


rule rmhost_report:
    input:
        reads = expand(os.path.join(config["results"]["rmhost"], "{{sample}}.rmhost{read}.fq.gz"),
                       read=[".1", ".2"] if IS_PE else "")
    output:
        os.path.join(config["results"]["report"]["rmhost"], "{sample}.rmhost.stats.tsv")
    params:
        fq_encoding = config["params"]["report"]["seqkit"]["fq_encoding"],
        sample_id = "{sample}"
    threads:
        config["params"]["report"]["seqkit"]["threads"]
    run:
        from metapi import reporter
        if IS_PE:
            shell("seqkit stats --all --basename --tabular \
                   --fq-encoding %s \
                   --out-file %s \
                   --threads %d %s" % (params.fq_encoding, output, threads, input))
            reporter.change(output[0], params.sample_id, "rmhost", "pe", ["fq1", "fq2"])
        else:
            shell("seqkit stats --all --basename --tabular \
                   --fq-encoding %s \
                   --out-file %s \
                   --threads %d %s" % (params.fq_encoding, output, threads, input))
            reporter.change(output[0], params.sample_id, "rmhost", "se", ["fq1"])


rule merge_rmhost_report:
    input:
        expand("{reportout}/{sample}.rmhost.stats.tsv",
               reportout=config["results"]["report"]["rmhost"],
               sample=_samples.index.unique())
    output:
        os.path.join(config["results"]["report"]["base_dir"], "rmhost.stats.tsv")
    run:
        from metapi import reporter
        reporter.merge(input, output[0], 8)
