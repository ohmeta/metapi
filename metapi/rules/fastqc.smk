def raw_reads(wildcards):
    if IS_PE:
        return [metapi.sampler.get_reads(SAMPLES, wildcards, "fq1"),
                metapi.sampler.get_reads(SAMPLES, wildcards, "fq2")]
    else:
        return [metapi.sampler.get_reads(SAMPLES, wildcards, "fq1")]


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
               sample=SAMPLES.index.unique())
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
            input_list.append(os.path.dirname(i))
        input_str = " ".join(input_list)
        shell("multiqc --outdir {params.outdir} --title fastqc --module fastqc %s 2> {log}" % input_str)


rule raw_report:
    input:
        unpack(raw_reads)
    output:
        os.path.join(config["results"]["report"]["raw"], "{sample}.raw.stats.tsv")
    params:
        fq_encoding = config["params"]["report"]["seqkit"]["fq_encoding"],
        sample_id = "{sample}"
    threads:
        config["params"]["report"]["seqkit"]["threads"]
    run:
        reads_num = len(input)
        if IS_PE:
            if reads_num == 2:
                shell("seqkit stats --all --basename --tabular \
                       --fq-encoding %s \
                       --out-file %s \
                       --threads %d %s" % (params.fq_encoding, output, threads, " ".join(input)))
                metapi.qcer.change(output[0], params.sample_id, "raw", "pe", ["fq1", "fq2"])
            else:
                r1_str = " ".join(input[0:reads_num//2])
                r2_str = " ".join(input[reads_num//2:])
                shell("cat %s | \
                       seqkit stats --all --basename --tabular \
                       --fq-encoding %s \
                       --out-file %s \
                       --threads %d" % (r1_str, params.fq_encoding, output[0] + ".1", threads))
                shell("cat %s | \
                       seqkit stats --all --basename --tabular \
                       --fq-encoding %s \
                       --out-file %s \
                       --threads %d" % (r2_str, params.fq_encoding, output[0] + ".2", threads))
                metapi.qcer.change(output[0] + ".1", params.sample_id, "raw", "pe", ["fq1"])
                metapi.qcer.change(output[0] + ".2", params.sample_id, "raw", "pe", ["fq2"])
                metapi.tooler.merge([output[0] + ".1", output[0] + ".2"], metapi.tooler.parse, 8, save=True, output=output[0])
                shell("rm -rf %s %s" % (output[0] + ".1", output[0] + ".2"))
        else:
            if reads_num == 1:
                shell("seqkit stats --all --basename --tabular \
                       --fq-encoding %s \
                       --out-file %s \
                       --threads %d %s" % (params.fq_encoding, output, threads))
                metapi.qcer.change(output[0], params.sample_id, "raw", "se", ["fq1"])
            else:
                r_str = " ".join(input)
                shell("cat %s | \
                       seqkit stats --all --basename --tabular \
                       --fq-encoding %s \
                       --out-file %s \
                       --threads %d" % (r_str, params.fq_encoding, output, threads))
                metapi.qcer.change(output[0], params.sample_id, "raw", "se", ["fq1"])


rule merge_raw_report:
    input:
        expand("{reportout}/{sample}.raw.stats.tsv",
               reportout=config["results"]["report"]["raw"],
               sample=SAMPLES.index.unique())
    output:
        os.path.join(config["results"]["report"]["base_dir"], "raw.stats.tsv")
    run:
        metapi.tooler.merge(input, metapi.tooler.parse, 8, save=True, output=output[0])
