rule build_asmfa_index:
    input:
        os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out/{sample}_{unit}.contigs.fa")
    output:
        expand("{assembly}/{{sample}}_{{unit}}.megahit_out/{{sample}}_{{unit}}.contigs.fa.{suffix}",
               assembly=config["results"]["assembly"],
               suffix=["amb", "ann", "bwt", "pac", "sa"])
    shell:
        "bwa index {input}"


rule align_reads_to_asmfa:
    input:
        reads = expand("{trim}/{{sample}}_{{unit}}.trimmed.{read}.fq.gz",
                       trim=config["results"]["trim"],
                       read=["1", "2"]),
        index = expand("{assembly}/{{sample}}_{{unit}}.megahit_out/{{sample}}_{{unit}}.contigs.fa.{suffix}",
                       assembly=config["results"]["assembly"],
                       suffix=["amb", "ann", "bwt", "pac", "sa"])
    output:
        flagstat = os.path.join(config["results"]["alignment"], "{sample}_{unit}.flagstat.txt"),
        bam = os.path.join(config["results"]["alignment"], "{sample}_{unit}.sorted.bam")
    params:
        bwa_mem_threads = config["params"]["alignment"]["bwa_mem_threads"],
        samtools_threads = config["params"]["alignment"]["samtools_threads"],
        prefix = os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out/{sample}_{unit}.contigs.fa")
    shell:
        "bwa mem -t {params.bwa_mem_threads} {params.prefix} {input.reads} | "
        "samtools view -@{params.samtools_threads} -hbS - | "
        "tee >(samtools flagstat -@{params.samtools_threads} - > {output.flagstat}) | "
        "samtools sort -@{params.samtools_threads} -o {output.bam} -"
