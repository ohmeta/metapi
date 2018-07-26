rule build_asmfa_index:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    output:
        expand("{assembly}/{{sample}}.megahit_out/{{sample}}.contigs.fa.gz.{suffix}",
               assembly=config["results"]["assembly"],
               suffix=["amb", "ann", "bwt", "pac", "sa"])
    shell:
        "bwa index {input}"

rule align_reads_to_asmfa:
    input:
        reads = assembly_inputs,
        index = expand("{assembly}/{{sample}}.megahit_out/{{sample}}.contigs.fa.gz.{suffix}",
                       assembly=config["results"]["assembly"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["alignment"], "{sample}.flagstat"),
        bam = os.path.join(config["results"]["alignment"], "{sample}.sorted.bam")
    params:
        bwa_mem_threads = config["params"]["alignment"]["bwa_mem_threads"],
        samtools_threads = config["params"]["alignment"]["samtools_threads"],
        prefix = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    shell:
        "bwa mem -t {params.bwa_mem_threads} {params.prefix} {input.reads} | "
        "samtools view -@{params.samtools_threads} -hbS - | "
        "tee >(samtools flagstat -@{params.samtools_threads} - > {output.flagstat}) | "
        "samtools sort -@{params.samtools_threads} -o {output.bam} -"
