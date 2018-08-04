def _get_fastq(wildcards, units, read_pair="fq1"):
    return units.loc[(wildcards.sample, wildcards.unit), [read_pair]].dropna()[0]

rule trimming_pe:
    input:
        fq1 = lambda wildcards: _get_fastq(wildcards, units, "fq1"),
        fq2 = lambda wildcards: _get_fastq(wildcards, units, "fq2")
    output:
        expand("{trim_dir}/{{sample}}_{{unit}}.trimmed.{read}.fq.gz",
               trim_dir=config["results"]["trim"],
               read=["1", "2", "single"])
    log:
        os.path.join(config["logs"]["trim"], "{sample}_{unit}.trimmed.log")
    params:
        qual_type = config["qual_type"],
    shell:
        "sickle pe -f {input.fq1} -r {input.fq2} "
        "-o {output[0]} -p {output[1]} -s {output[2]} "
        "--gzip-output -t {params.qual_type} {log}"


rule trimming_se:
    input:
        fq = lambda wildcards: get_fastq(wildcards, units)
    output:
        fq = os.path.join(config["results"]["trim"], "{sample}}_{unit}.trimmed.fq.gz")
    log:
        os.path.join(config["logs"]["trim"], "{sample}_{unit}.trimmed.log")
    params:
        qual_type = config["qual_type"],
    shell:
        "sickle se -f {input.fq} -o {output.fq} "
        "--gzip-output -t {params.qual_type} {log}"
