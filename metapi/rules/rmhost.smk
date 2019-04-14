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
            r1 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz"),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                           suffix=["amb", "ann", "bwt", "pac", "sa"])
        output:
            flagstat = protected(os.path.join(config["results"]["rmhost"], "{sample}.rmhost.flagstat.txt")),
            r1 = protected(os.path.join(config["results"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
            r2 = protected(os.path.join(config["results"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
        log:
            os.path.join(config["logs"]["rmhost"], "{sample}.bwa.rmhost.log")
        params:
            minimum_seed_length = config["params"]["rmhost"]["bwa"]["minimum_seed_length"],
            save_bam = "true" if config["params"]["rmhost"]["bwa"]["save_bam"] else "false",
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bam = os.path.join(config["results"]["rmhost"], "{sample}.bwa.host.sorted.bam")
        threads:
            config["params"]["rmhost"]["bwa"]["threads"]
        shell:
            '''
            if {params.save_bam}; then
                bwa mem -k {params.minimum_seed_length} -t {threads} {params.index_prefix} {input.r1} {input.r2} | \
                tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                tee >(samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.r1} -2 {output.r2} -) | \
                samtools sort -@{threads} -O BAM -o {params.bam} - 2>{log}
            else
                bwa mem -k {params.minimum_seed_length} -t {threads} {params.prefix} {input.r1} {input.r2} | \
                tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.r1} -2 {output.r2} - 2>{log}
            fi
            '''

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
            r1 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz"),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            flagstat = protected(os.path.join(config["results"]["rmhost"], "{sample}.rmhost.flagstat.txt")),
            r1 = protected(os.path.join(config["results"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
            r2 = protected(os.path.join(config["results"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
        log:
            os.path.join(config["logs"]["rmhost"], "{sample}.bowtie2.rmhost.log")
        params:
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"],
            additional_params = config["params"]["rmhost"]["bowtie2"]["additional_params"],
            bam = os.path.join(config["results"]["rmhost"], "{sample}.bowtie2.host.sorted.bam"),
            save_bam = "true" if config["params"]["rmhost"]["bwa"]["save_bam"] else "false",
        threads:
            config["params"]["rmhost"]["bowtie2"]["threads"]
        shell:
            '''
            if {params.save_bam}; then
                bowtie2 --threads {threads} -x {params.index_prefix} \
                -1 {input.r1} -2 {input.r2} {params.additional_params} 2> {log} | \
                tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                tee >(samtools sort -@{threads} -O BAM -o {params.bam}) | \
                samtools view -@{threads} -SF4 - | awk -F'[/\t]' '{{print $1}}' | sort | uniq | \
                tee >(awk '{{print $0 "/1"}}' - | seqtk subseq -r {input.r1} - | pigz -p {threads} -c > {output.r1}) | \
                awk '{{print $0 "/2"}}' - | seqtk subseq -r {input.r2} - | pigz -p {threads} -c > {output.r2}
            else
                bowtie2 --threads {threads} -x {params.index_prefix} \
                -1 {input.r1} -2 {input.r2} {params.additional_params} 2> {log} | \
                tee >(samtools flagstat -@{threads} - > {output.flagstat}) | \
                samtools view -@{threads} -SF4 - | awk -F'[/\t]' '{{print $1}}' | sort | uniq | \
                tee >(awk '{{print $0 "/1"}}' - | seqtk subseq -r {input.r1} - | pigz -p {threads} -c > {output.r1}) | \
                awk '{{print $0 "/2"}}' - | seqtk subseq -r {input.r2} - | pigz -p {threads} -c > {output.r2}
            fi
            '''
