rule trimming_pe:
    input:
        fq1 = lambda wildcards: samples[wildcards.sample][0],
        fq2 = lambda wildcards: samples[wildcards.sample][1]
    output:
        expand("{trim}/{sample}.trimmed.{read}.fq.gz",
               trim=config["results"]["trim"],
               sample=samples.keys(),
               read=["1", "2", "single"])
    log:
        os.path.join(config["logs"]["trim"], "{sample}.trimmed.log")
    params:
        qual_type = config["params"]["trim"]["qual_type"],
        qual_cutoff = config["params"]["trim"]["qual_cutoff"],
        length_cutoff = config["params"]["trim"]["length_cutoff"]
    shell:
        "sickle pe -f {input.fq1} -r {input.fq2} "
        "-o {output[0]} -p {output[1]} -s {output[2]} "
        "--gzip-output -t {params.qual_type} -q {params.qual_cutoff} -l {params.length_cutoff} {log}"