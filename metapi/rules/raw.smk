def raw_short_reads(wildcards):
    if READS_FORMAT == "fastq":
        if config["params"]["simulate"]["do"]:
            return [metapi.get_reads(SAMPLES, wildcards, "fq1")[0],
                    metapi.get_reads(SAMPLES, wildcards, "fq2")[0]]
        else:
            if IS_PE:
                if IS_INTERLEAVED:
                    return [metapi.get_reads(SAMPLES, wildcards, "fq1")]
                else:
                    return [metapi.get_reads(SAMPLES, wildcards, "fq1"),
                            metapi.get_reads(SAMPLES, wildcards, "fq2")]
            else:
                return [metapi.get_reads(SAMPLES, wildcards, "fq1")]
    elif READS_FORMAT == "sra":
        return [metapi.get_reads(SAMPLES, wildcards, "sra")]


def raw_long_reads(wildcards):
    if READS_FORMAT == "fastq":
        return [metapi.get_reads(SAMPLES, wildcards, "fq_long")]
    else:
        print("Don't support SRA format now, exit.")
        sys.exit(1)


def reads_suffix():
    if HAVE_LONG:
        return [".1", ".2", ".long"]
    elif IS_PE:
        return [".1", ".2"]
    else:
        return [""]


def short_reads_suffix():
    if IS_PE:
        return [".1", ".2"]
    else:
        return [""]

 
def long_reads_suffix():
    return [".long"]


if config["params"]["raw"]["do"]:
    rule prepare_short_reads:
        input:
            unpack(raw_short_reads)
        output:
            reads = expand(
                os.path.join(
                    config["output"]["raw"],
                    "short_reads/{{sample}}/{{sample}}.raw{read}.fq.gz"),
                read=short_reads_suffix()) \
                if config["params"]["raw"]["save_reads"] else \
                temp(expand(
                    os.path.join(
                        config["output"]["raw"],
                        "short_reads/{{sample}}/{{sample}}.raw{read}.fq.gz"),
                    read=short_reads_suffix()))
        params:
            output_dir = os.path.join(config["output"]["raw"], "short_reads/{sample}"),
            is_pe = IS_PE,
            is_interleaved = IS_INTERLEAVED,
            reads_format = READS_FORMAT,
            check_paired = config["params"]["raw"]["check_paired"]
        threads:
            config["params"]["raw"]["threads"]
        conda:
            config["envs"]["prepare"]
        log:
            os.path.join(config["output"]["raw"], "logs/{sample}_prepare.log")
        script:
            "../wrappers/preprocess_raw.py"


    rule prepare_short_reads_all:
        input:
            expand(os.path.join(
                config["output"]["raw"],
                "short_reads/{sample}/{sample}.raw{read}.fq.gz"),
                read=short_reads_suffix(),
                sample=SAMPLES_ID_LIST)

else:
    rule prepare_short_reads_all:
        input:


if config["params"]["raw"]["do"]:
    if HAVE_LONG:
        rule prepare_long_reads:
            input:
                unpack(raw_long_reads)
            output:
                reads = expand(
                    os.path.join(
                        config["output"]["raw"],
                        "long_reads/{{sample}}/{{sample}}.raw{read}.fq"),
                    read=long_reads_suffix()) \
                    if config["params"]["raw"]["save_reads"] else \
                    temp(expand(
                        os.path.join(
                            config["output"]["raw"],
                            "long_reads/{{sample}}/{{sample}}.raw{read}.fq"),
                        read=long_reads_suffix()))
            run:
                reads_num = len(input)

                if READS_FORMAT == "fastq":
                    if reads_num == 1:
                        if input[0].endswith(".gz"):
                            shell('''gzip -dc %s > %s''' % (input[0], output.reads[0]))
                        else:
                            os.symlink(os.path.realpath(input[0]), output.reads[0])
                    else:
                        for i in input:
                            if i.endswith(".gz"):
                                shell('''gzip -dc %s >> %s''' % (i, output.reads[0]))
                            else:
                                shell('''cat %s >> %s''' % (i, output.reads[0]))


        rule prepare_long_reads_all:
            input:
                expand(os.path.join(
                    config["output"]["raw"],
                    "long_reads/{sample}/{sample}.raw{read}.fq"),
                    read=long_reads_suffix(),
                    sample=SAMPLES_ID_LIST)

    else:
        rule prepare_long_reads_all:
            input:

else:
    rule prepare_long_reads_all:
        input:


rule prepare_reads_all:
    input:
        rules.prepare_short_reads_all.input,
        rules.prepare_long_reads_all.input


def get_reads(wildcards, step, have_single=False, have_long=False):
    read = short_reads_suffix()
    if IS_PE and have_single:
        read += [".single"]

    short_reads = expand(os.path.join(
        config["output"][step],
        "short_reads/{sample}/{sample}.{step}{read}.fq.gz"),
                         step=step,
                         read=read,
                         sample=wildcards.sample)

    long_reads = expand(os.path.join(
        config["output"]["raw"],
        "long_reads/{sample}/{sample}.{step}{read}.fq"),
                        step="raw",
                        read=long_reads_suffix(),
                        sample=wildcards.sample)

    if have_long:
        return short_reads + long_reads
    else:
        return short_reads


def get_reads_(wildcards, step, have_single=False):
    read = ""
    if IS_PE:
        read = short_reads_suffix()
        if have_single:
            read += [".single"]

    return expand(
        os.path.join(
            config["output"][step],
            "short_reads/{sample_}/{sample_}.{step}{read}.fq.gz"),
        step=step,
        read=read,
        sample_=wildcards.sample_)


def get_short_reads_list(step, samples_id_list):
    if IS_PE:
        return [
            expand(
                os.path.join(
                    config["output"][step],
                    "short_reads/{sample}/{sample}.{step}.1.fq.gz"),
                step=step,
                sample=samples_id_list),
            expand(
                os.path.join(
                    config["output"][step],
                    "short_reads/{sample}/{sample}.{step}.2.fq.gz"),
                step=step,
                sample=samples_id_list)]
    else:
        return [
            expand(
                os.path.join(
                    config["output"][step],
                    "short_reads/{sample}/{sample}.{step}.fq.gz"),
                step=step,
                sample=samples_id_list)]


def get_long_reads_list(step, samples_id_list):
    return [
        expand(
            os.path.join(
                config["output"][step],
                "long_reads/{sample}/{sample}.{step}.long.fq.gz"),
            step=step,
            sample=samples_id_list)]


if config["params"]["raw"]["fastqc"]["do"]:
    rule raw_fastqc:
        input:
            lambda wildcards: get_reads(wildcards, "raw", False)
        output:
            directory(os.path.join(
                config["output"]["raw"],
                "fastqc/{sample}.fastqc.out"))
        conda:
            config["envs"]["fastqc"]
        threads:
            config["params"]["raw"]["threads"]
        log:
            os.path.join(config["output"]["raw"], "logs/{sample}.fastqc.log")
        shell:
            '''
            mkdir -p {output}
            fastqc \
            --outdir {output} \
            --threads {threads} \
            --format fastq \
            {input} \
            2> {log}
            '''


    rule raw_fastqc_multiqc:
        input:
            expand(os.path.join(
                config["output"]["raw"],
                "fastqc/{sample}.fastqc.out"),
                   sample=SAMPLES_ID_LIST)
        output:
            html = os.path.join(
                config["output"]["raw"],
                "report/fastqc_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["raw"],
                "report/fastqc_multiqc_report_data"))
        conda:
            config["envs"]["multiqc"]
        params:
            output_dir = os.path.join(config["output"]["raw"],
                                      "report")
        log:
            os.path.join(config["output"]["raw"], "logs/multiqc_fastqc.log")
        threads: 
            1
        shell:
            '''
            multiqc \
            --outdir {params.output_dir} \
            --title fastqc \
            --module fastqc {input} \
            2> {log}
            '''


    rule raw_fastqc_all:
        input:
            expand([
                os.path.join(
                    config["output"]["raw"],
                    "fastqc/{sample}.fastqc.out"),
                os.path.join(
                    config["output"]["raw"],
                    "report/fastqc_multiqc_report.html"),
                os.path.join(
                    config["output"]["raw"],
                    "report/fastqc_multiqc_report_data")],
                   sample=SAMPLES_ID_LIST)

else:
    rule raw_fastqc_all:
        input:


if config["params"]["qcreport"]["do"]:
    rule raw_report:
        input:
            lambda wildcards: get_reads(wildcards, "raw")
        output:
            temp(os.path.join(config["output"]["raw"],
                              "report/stats/{sample}_raw_stats.tsv.raw"))
        conda:
            config["envs"]["report"]
        params:
            fq_encoding = config["params"]["fq_encoding"]
        log:
            os.path.join(config["output"]["raw"],
                         "logs/{sample}.seqkit.log")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        shell:
            '''
            seqkit stats \
            --all \
            --basename \
            --tabular \
            --fq-encoding {params.fq_encoding} \
            --out-file {output} \
            --threads {threads} \
            {input} 2> {log}
            '''


    rule raw_report_refine:
        input:
            os.path.join(config["output"]["raw"],
                         "report/stats/{sample}_raw_stats.tsv.raw")
        output:
            os.path.join(config["output"]["raw"],
                         "report/stats/{sample}_raw_stats.tsv.gz")
        params:
            sample_id = "{sample}"
        threads:
            1
        run:
            if IS_PE:
                metapi.change(str(input), str(output), params.sample_id, "raw",
                              "pe", ["fq1", "fq2"])
            else:
                metapi.change(str(input), str(output), params.sample_id, "raw",
                              "se", ["fq1"])


    rule raw_report_merge:
        input:
            expand(
                os.path.join(config["output"]["raw"],
                             "report/stats/{sample}_raw_stats.tsv.gz"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "raw_stats.tsv.gz")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            metapi.merge(input, metapi.parse, threads, output=output[0])


    rule raw_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "raw_stats.tsv.gz")

else:
    rule raw_report_all:
        input:


rule raw_all:
    input:
        #rules.prepare_reads_all.input,
        rules.raw_fastqc_all.input,
        rules.raw_report_all.input


localrules:
    prepare_short_reads_all,
    prepare_long_reads_all,
    raw_fastqc_all,
    raw_report_all,
    raw_all