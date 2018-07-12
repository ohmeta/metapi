rule fastqc:
    input:
        lambda wildcards: samples[wildcards.sample]
    output:
        outdir = config["results"]["raw"]["fastqc"],
        outfile = expand("{fastqc}/{sample}_{read}_fastqc.{out}",
                         fastqc=config["results"]["raw"]["fastqc"],
                         sample=samples.keys(),
                         read=["1", "2"],
                         out=["html", "zip"])
    shell:
        "fastqc -o {output.outdir} -f fastq {input}"