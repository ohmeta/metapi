rule fastqc:
    input:
        r1 = lambda wildcards: _get_raw_fastq(wildcards, "fq1"),
        r2 = lambda wildcards: _get_raw_fastq(wildcards, "fq2")
    output:
        outfile = expand("{fastqc}/{{sample}}_{read}_fastqc.{out}",
                         fastqc=config["results"]["raw"]["fastqc"],
                         read=["1", "2"],
                         out=["html", "zip"])
    params:
        outdir = config["results"]["raw"]["fastqc"]
    shell:
        "fastqc -o {params.outdir} -f fastq {input.r1} {input.r2}"
