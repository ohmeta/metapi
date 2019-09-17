def raw_reads(wildcards):
    if IS_PE:
        return [sample.get_reads(_samples, wildcards, "fq1"),
                sample.get_reads(_samples, wildcards, "fq2")]
    else:
        return [sample.get_reads(_samples, wildcards, "fq1")]


rule fastqc:
    input:
        unpack(raw_reads)
    output:
        done = os.path.join(config["results"]["raw"]["fastqc"], "{sample}/done")
    params:
        outdir = os.path.join(config["results"]["raw"]["fastqc"], "{sample}")
    threads:
        config["params"]["fastqc"]["threads"]
    log:
        os.path.join(config["logs"]["raw"]["fastqc"], "{sample}_fastqc.log")
    shell:
        '''
        fastqc -o {params.outdir} -t {threads} -f fastq {input} 2> {log}
        echo done > {output.done}
        '''

rule multiqc_fastqc:
    input:
        expand("{fastqc}/{sample}/done",
               fastqc=config["results"]["raw"]["fastqc"],
               sample=_samples.index.unique())
    output:
        html = os.path.join(config["results"]["raw"]["multiqc"], "fastqc_multiqc_report.html"),
        data_dir = directory(os.path.join(config["results"]["raw"]["multiqc"], "fastqc_multiqc_report_data"))
    params:
        outdir = config["results"]["raw"]["multiqc"]
    log:
        os.path.join(config["logs"]["raw"]["multiqc"], "multiqc_fastqc.log")
    run:
        input_list = []
        for i in input:
            input_list.append(os.path.basename(i))
        input_str = " ".join(input_list)
        shell("multiqc --outdir {params.outdir} --title fastqc --module fastqc %s 2> {log}" % input_str)
