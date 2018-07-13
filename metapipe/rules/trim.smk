rule trim_pe:
    input:
        r1 = lambda wildcards: get_fastq(wildcards)[0],
        r2 = lambda wildcards: get_fastq(wildcards)[1]
    output:
        expand("{trim}/{{sample}}.trimmed.{read}.fq.gz",
               trim=config["results"]["trim"],
               read=["1", "2", "single"])
    params:
        qual_type = config["params"]["trim"]["qual_type"],
        qual_cutoff = config["params"]["trim"]["qual_cutoff"],
        length_cutoff = config["params"]["trim"]["length_cutoff"]
    shell:
        "sickle pe -f {input.r1} -r {input.r2} "
        "-o {output[0]} -p {output[1]} -s {output[2]} "
        "--gzip-output -t {params.qual_type} -q {params.qual_cutoff} -l {params.length_cutoff}"