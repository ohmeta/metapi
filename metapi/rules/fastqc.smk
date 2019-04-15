rule fastqc:
    input:
        r1 = lambda wildcards: sample.get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: sample.get_sample_id(_samples, wildcards, "fq2")
    output:
        outfile = expand("{fastqc}/{{sample}}_{read}_fastqc.{out}",
                         fastqc=config["results"]["raw"]["fastqc"],
                         read=["1", "2"],
                         out=["html", "zip"])
    params:
        outdir = config["results"]["raw"]["fastqc"]
    log:
        os.path.join(config["logs"]["raw"]["fastqc"], "{sample}_fastqc.log")
    shell:
        "fastqc -o {params.outdir} -f fastq {input.r1} {input.r2} 2> {log}"

rule multiqc_fastqc:
    input:
        expand("{fastqc}/{sample}_{read}_fastqc.zip",
               fastqc=config["results"]["raw"]["fastqc"],
               sample=_samples.index,
               read=["1", "2"])
    output:
        html = os.path.join(config["results"]["raw"]["multiqc"], "fastqc_multiqc_report.html"),
        data_dir = directory(os.path.join(config["results"]["raw"]["multiqc"], "fastqc_multiqc_report_data"))
    params:
        outdir = config["results"]["raw"]["multiqc"]
    log:
        os.path.join(config["logs"]["raw"]["multiqc"], "multiqc_fastqc.log")
    shell:
        '''
        multiqc --outdir {params.outdir}  --title fastqc --module fastqc {input} 2> {log}
        '''
