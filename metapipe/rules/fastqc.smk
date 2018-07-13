def get_fastq(wildcards):
    return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

rule fastqc:
    input:
        get_fastq
    output:
        outfile = expand("{fastqc}/{{sample}}_{read}_fastqc.{out}",
                         fastqc=config["results"]["raw"]["fastqc"],
                         read=["1", "2"],
                         out=["html", "zip"])
    params:
        outdir = config["results"]["raw"]["fastqc"]
    shell:
        """
        echo {input}
        fastqc -o {params.outdir} -f fastq {input}
        """