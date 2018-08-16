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
        reads = expand("{trimming}/{{sample}}.trimmed.{read}.fq.gz",
                       trimming=config["results"]["trimming"],
                       read=["1", "2"]),
        index = expand("{prefix}.{suffix}",
                       prefix=config["results"]["host"]["prefix"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["rmhost"], "{sample}.flagstat.txt"),
        reads = expand("{rmhost}/{{sample}}.rmhost.{read}.fq.gz",
                       rmhost=config["results"]["rmhost"],
                       read=["1", "2"])
    params:
        threads = config["params"]["rmhost"]["threads"],
        prefix = config["results"]["host"]["prefix"]
    shell:
        "bwa mem -t {threads} {params.prefix} {input.reads} | "
        "tee >(samtools flagstat -@{threads} - > {output.flagstat}) | "
        "samtools fastq -@{threads} -f 12 -n -1 {output.reads[0]} -2 {output.reads[1]} -"
