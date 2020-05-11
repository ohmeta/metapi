rule prepare_reads:
    input:
        unpack(raw_reads)
    output:
        expand(
            os.path.join(config["output"]["raw"],
                         "short_reads/{{sample}}/{{sample}}.raw{read}.fq.gz"),
            read=[".1", ".2"] if IS_PE else "")
    params:
        output_dir = os.path.join(config["output"]["raw"],
                                  "short_reads/{sample}")
    threads:
        config["params"]["raw"]["threads"]
    run:
        reads_num = len(input)

        if READS_FORMAT == "fastq":
            if IS_PE:
                if not IS_INTERLEAVED:
                    if reads_num == 2:
                        os.symlink(os.path.realpath(input[0]), output[0])
                        os.symlink(os.path.realpath(input[1]), output[1])
                    else:
                        shell('''cat %s > %s''' % (" ".join(input[0:reads_num//2]), output[0]))
                        shell('''cat %s > %s''' % (" ".join(input[reads_num//2:]), output[1]))
                else:
                    shell(
                        '''
                        cat {input} | \
                        tee >(seqtk seq -1 - | pigz -c -p {threads} > {output[0]}) | \
                        seqtk seq -2 - | pigz -c -p {threads} > {output[1]}
                        ''')
            else:
                if reads_num == 1:
                    os.symlink(os.path.realpath(input[0]), output[0])
                else:
                    shell('''cat {input} > {output}''')

        elif READS_FORMAT == "sra":
            reads_direction = str("+"),
            header_format = str("@$ac-$si/$ri")

            if reads_num == 1:
                sra_id = os.path.basename(input[0]).split(".")[0]
                shell(
                    '''
                    fastq-dump \
                    --gzip \
                    --split-3 \
                    --defline-qual '%s' \
                    --defline-seq '%s' \
                    --outdir {params.output_dir} \
                    {input[0]}
                    ''' % (reads_direction, header_format))

                shell('''mv {params.output_dir}/%s_1.fastq.gz {output[0]}''' % sra_id)
                shell('''mv {params.output_dir}/%s_2.fastq.gz {output[1]}''' % sra_id)
                shell('''rm -rf {params.output_dir}/%s.fastq.gz''' % sra_id)
            else:
                r1_list = []
                r2_list = []
                for sra_file in input:
                    sra_id = os.path.basename(sra_file).split(".")[0]
                    r1_list.append(os.path.join(params.output_dir,
                                                sra_id + "_1.fastq.gz"))
                    r2_list.append(os.path.join(params.output_dir,
                                                sra_id + "_2.fastq.gz"))
                    shell(
                        '''
                        fastq-dump \
                        --gzip \
                        --split-3 \
                        --defline-qual '%s' \
                        --defline-seq '%s' \
                        --outdir {params.output_dir} %s
                        ''' % reads_direction, header_format, sra_file)

                    shell('''rm -rf {params.output_dir}/%s.fastq.gz''' % sra_id)

                    r1_str = " ".join(r1_list)
                    r2_str = " ".join(r2_list)
                    shell('''cat %s > %s''' % (r1_str, output[0]))
                    shell('''cat %s > %s''' % (r2_str, output[1]))
                    shell('''rm -rf %s''' % r1_str)
                    shell('''rm -rf %s''' % r2_str)


rule prepare_reads_all:
    input:
         expand(
             os.path.join(config["output"]["raw"],
                          "short_reads/{sample}/{sample}.raw{read}.fq.gz"),
             read=[".1", ".2"] if IS_PE else "",
             sample=SAMPLES.index.unique())

        
def get_reads(wildcards, step, have_single=False):
    read = ""
    if IS_PE:
        read = [".1", ".2"]
        if have_single:
            read += [".single"]

    return expand(
        os.path.join(
            config["output"][step],
            "short_reads/{sample}/{sample}.{step}{read}.fq.gz"),
        step=step,
        read=read,
        sample=wildcards.sample)


def get_reads_(wildcards, step, have_single=False):
    read = ""
    if IS_PE:
        read = [".1", ".2"]
        if have_single:
            read += [".single"]

    return expand(
        os.path.join(
            config["output"][step],
            "short_reads/{sample_}/{sample_}.{step}{read}.fq.gz"),
        step=step,
        read=read,
        sample_=wildcards.sample_)


def get_reads_list(step):
    if IS_PE:
        return [
            expand(
                os.path.join(
                    config["output"][step],
                    "short_reads/{sample}/{sample}.{step}.1.fq.gz"),
                step=step,
                sample=SAMPLES.index.unique()),
            expand(
                os.path.join(
                    config["output"][step],
                    "short_reads/{sample}/{sample}.{step}.2.fq.gz"),
                step=step,
                sample=SAMPLES.index.unique())]
    else:
        return [
            expand(
                os.path.join(
                    config["output"][step],
                    "short_reads/{sample}/{sample}.{step}.fq.gz"),
                step=step,
                sample=SAMPLES.index.unique())]


if config["params"]["raw"]["fastqc"]["do"]:
    rule raw_fastqc:
        input:
            lambda wildcards: get_reads(wildcards, "raw", False)
        output:
            directory(os.path.join(
                config["output"]["raw"],
                "fastqc/{sample}.fastqc.out"))
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
                   sample=SAMPLES.index.unique())
        output:
            html = os.path.join(
                config["output"]["raw"],
                "report_multiqc/fastqc_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["raw"],
                "report_multiqc/fastqc_multiqc_report_data"))
        params:
            output_dir = os.path.join(config["output"]["raw"],
                                      "report_multiqc")
        log:
            os.path.join(config["output"]["raw"], "logs/multiqc_fastqc.log")
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
                    "report_multiqc/fastqc_multiqc_report.html"),
                os.path.join(
                    config["output"]["raw"],
                    "report_multiqc/fastqc_multiqc_report_data")],
                   sample=SAMPLES.index.unique())

else:
    rule raw_fastqc_all:
        input:


rule raw_all:
    input:
        rules.prepare_reads_all.input,
        rules.raw_fastqc_all.input,
