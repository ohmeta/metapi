rule trimming:
    input:
        fastq1 = lambda wildcards: config["R1_raw"][wildcards.sample],
        fastq2 = lambda wildcards: config["R2_raw"][wildcards.sample]
    output:
        expand("{filter_dir}/{{sample}}.clean.{ext}.", filter_dir=config["filter_reads_dir"], ext=["1.fq.gz", "2.fq.gz", "single.fq.gz", "stat_out"])
    params:
        oa_filter = config["oa_filter"],
        phread_quality_system = config["phread_quality_system"],
        oa_min_length = config["oa_min_length"],
        out_dir = config["filter_reads_dir"],
        seed_oa = config["seed_oa"],
        frag_oa = config["frag_oa"]
    run:
        prefix = os.path.join("{params.out_dir}", "{sample}")
        shell("""perl {params.oa_filter} {input.R1},{input.R2} {prefix} {params.phread_quality_system} {params.oa_min_length} {params.seed_oa} {params.frag_oa}""")
