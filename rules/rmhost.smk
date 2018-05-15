rule build_host_index:
    input:
        config["host_index"]["fasta"]
    output:
        expand("{prefix}.{suffix}",
               prefix=config["host_index"]["prefix"],
               suffix=["amb", "ann", "bwt", "pac", "sa"])
    params:
        prefix = config["host_index"]["prefix"]
    shell:
        "bwa index {input} -p {params.prefix}"


rule rmhost_pe_pair:
    input:
        reads = expand("{trim_dir}/{{sample}}_{{unit}}.trimmed.{read}.fq.gz",
                       trim_dir=config["results"]["trim"],
                       read=["1", "2"]),
        index = expand("{prefix}.{suffix}",
                       prefix=config["host_index"]["prefix"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["rmhost"], "{sample}_{unit}.flagstat.txt"),
        reads = expand("{rmhost_dir}/{{sample}}_{{unit}}.rmhost.{read}.fq.gz",
                       rmhost_dir=config["results"]["rmhost"],
                       read=["1", "2"])
    params:
        bwa_mem_threads = config["params"]["rmhost"]["bwa_mem_threads"],
        samtools_threads = config["params"]["rmhost"]["samtools_threads"],
        prefix = config["host_index"]["prefix"]
    shell:
        "bwa mem -t {params.bwa_mem_threads} {params.prefix} {input.reads} | "
        "tee >(samtools flagstat -@{params.samtools_threads} - > {output.flagstat}) | "
        "samtools fastq -@{params.samtools_threads} -f 12 -n -1 {output.reads[0]} -2 {output.reads[1]} -"
