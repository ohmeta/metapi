rule rmhost__aa:
    input:
        r1 = os.path.join(config["results"]["trimming"], "trimmed.{sample}_1.fq.gz"),
        r2 = os.path.join(config["results"]["trimming"], "trimmed.{sample}_2.fq.gz")
    output:
        flagstat = os.path.join(config["results"]["rmhost"], "{sample}_rmhostaa_flagstat.txt"),
        r1 = os.path.join(config["results"]["rmhost"], "rmhostaa.{sample}_1.fq.gz"),
        r2 = os.path.join(config["results"]["rmhost"], "rmhostaa.{sample}_2.fq.gz"),
        bam = os.path.join(config["results"]["rmhost"], "{sample}_host_sorted.bam")
    log:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhostaa.log")
    params:
        save_bam = config["params"]["rmhost"]["save_bam"],
        prefix = config["results"]["host"]["prefix"]
    threads:
        config["params"]["rmhost"]["threads"]
    shell:
        '''
        bwa mem -t {threads} {params.prefix} {input.r1} {input.r2} |
        tee >(samtools flagstat -@{threads} - > {output.flagstat}) |
        tee >(samtools fastq -@{threads} -N -f 12 -F 256 -1 {output.r1} -2 {output.r2} -) |
        samtools sort -@{threads} -o {output.bam} - 2>{log}
        '''
