if not os.path.exists(os.path.join(config["results"]["host"]["prefix"], "bwt")):
    rule build_host_index:
        input:
            config["results"]["host"]["fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["results"]["host"]["prefix"],
                   suffix=["amb", "ann", "bwt", "pac", "sa"])
        params:
            prefix = config["results"]["host"]["prefix"]
        shell:
            "bwa index {input} -p {params.prefix}"


rule rmhost:
    input:
        r1 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz"),
        r2 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz")
    output:
        flagstat = os.path.join(config["results"]["rmhost"], "{sample}.flagstat.txt"),
        r1 = os.path.join(config["results"]["rmhost"], "rmhost.{sample}.rmhost.1.fq.gz"),
        r2 = os.path.join(config["results"]["rmhost"], "rmhost.{sample}.rmhost.2.fq.gz"),
        bam = os.path.join(config["results"]["rmhost"], "{sample}.host.sorted.bam")
    log:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.log")
    params:
        save_bam = config["params"]["rmhost"]["save_bam"],
        prefix = config["results"]["host"]["prefix"]
    threads:
        config["params"]["rmhost"]["threads"]
    run:
        if {params.save_bam}:
            shell(
                '''
                bwa mem -t {threads} {params.prefix} {input.r1} {input.r2} |
                tee >(samtools flagstat -@{threads} - > {output.flagstat}) |
                tee >(samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.r1} -2 {output.r2} -) |
                samtools sort -@{threads} -o {output.bam} - 2>{log}
                ''')
        else:
            shell(
                '''
                bwa mem -t {threads} {params.prefix} {input.r1} {input.r2} |
                tee >(samtools flagstat -@{threads} - > {output.flagstat}) |
                samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.r1} -2 {output.r2} - 2>{log}
                ''')

