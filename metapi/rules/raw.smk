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
            output_dir = os.path.join(config["output"]["raw"],
                                   "short_reads/{sample}"),
            interleaved = config["params"]["interleaved"]
        threads:
            config["params"]["raw"]["threads"]
        log:
            os.path.join(config["output"]["raw"], "logs/{sample}_prepare.log")
        run:
            import filecmp
            import gzip

            reads_num = len(input)

            if READS_FORMAT == "fastq":
                if IS_PE:
                    if not params.interleaved:
                        if reads_num == 2:
                            os.symlink(os.path.realpath(input[0]), f'''{output.reads[0]}.temp.gz''')
                            os.symlink(os.path.realpath(input[1]), f'''{output.reads[1]}.temp.gz''')
                        else:
                            shell(f'''cat {" ".join(input[0:reads_num//2])} > {output.reads[0]}.temp.gz 2>> {log}''')
                            shell(f'''cat {" ".join(input[reads_num//2:])}  > {output.reads[1]}.temp.gz 2>> {log}''')

                        shell(f'''seqkit seq -ni {output.reads[0]}.temp.gz | sed 's#/1$##g' > {params.output_dir}/id.list.1 2>> {log}''') 
                        shell(f'''seqkit seq -ni {output.reads[1]}.temp.gz | sed 's#/2$##g' > {params.output_dir}/id.list.2 2>> {log}''') 

                        if filecmp.cmp(f'''{params.output_dir}/id.list.1''', f'''{params.output_dir}/id.list.2'''): 
                            shell(f'''mv {output.reads[0]}.temp.gz {output.reads[0]} 2>> {log}''')
                            shell(f'''mv {output.reads[1]}.temp.gz {output.reads[1]} 2>> {log}''')
                        else:
                            shell(
                                '''
                                sort -T {params.output_dir} {params.output_dir}/id.list.1 {params.output_dir}/id.list.2 | \
                                uniq -c | \
                                awk '$1==2{{print $2}}' > {params.output_dir}/id.list.paired 2>> {log}
                                ''')

                            oneline = gzip.open(f'''{output.reads[0]}.temp.gz''', 'rt').readline().strip().split()[0]
                            if "/1" in oneline:
                                shell(
                                    '''
                                    seqkit grep -f <(awk '{{print $0 "/1"}}' {params.output_dir}/id.list.paired) {output.reads[0]}.temp.gz -o {output.reads[0]} 2>> {log}
                                    seqkit grep -f <(awk '{{print $0 "/2"}}' {params.output_dir}/id.list.paired) {output.reads[1]}.temp.gz -o {output.reads[1]} 2>> {log}
                                    ''')
                            else:
                                shell(
                                    f'''
                                    seqkit grep -f {params.output_dir}/id.list.paired {output.reads[0]}.temp.gz -o {output.reads[0]} 2>> {log}
                                    seqkit grep -f {params.output_dir}/id.list.paired {output.reads[1]}.temp.gz -o {output.reads[1]} 2>> {log}
                                    ''')

                            shell(f'''rm -rf {output.reads[0]}.temp.gz 2>> {log}''')
                            shell(f'''rm -rf {output.reads[1]}.temp.gz 2>> {log}''')
                            shell(f'''rm -rf {params.output_dir}/id.list.paired 2>> {log}''')

                        shell(f'''rm -rf {params.output_dir}/id.list.1 2>> {log}''')
                        shell(f'''rm -rf {params.output_dir}/id.list.2 2>> {log}''')

                    else:
                        shell(
                            '''
                            cat {input} | \
                            tee >(seqtk seq -1 - | pigz -c -p {threads} > {output.reads[0]}) | \
                            seqtk seq -2 - | pigz -c -p {threads} > {output.reads[1]} 2>> {log}
                            ''')
                else:
                    if reads_num == 1:
                        os.symlink(os.path.realpath(input[0]), output.reads[0])
                    else:
                        shell('''cat {input} > {output.reads[0]} 2>> {log}''')

            elif READS_FORMAT == "sra":
                if reads_num == 1:
                    sra_file = os.path.basename(input[0])
                    shell(
                        f'''
                        rm -rf {params.output_dir}/{sra_file}* 2>> {log}
                        rm -rf {params.output_dir}.{sra_file}.temp 2>> {log}

                        fasterq-dump \
                        --threads {threads} \
                        --split-3 \
                        --temp {params.output_dir}.{sra_file}.temp \
                        --outdir {params.output_dir} {input[0]} 2>>{log}

                        rm -rf {params.output_dir}.{sra_file}.temp 2>> {log}
                        pigz --processes {threads} {params.output_dir}/{sra_file}_1.fastq 2>> {log}
                        pigz --processes {threads} {params.output_dir}/{sra_file}_2.fastq 2>> {log}
                        rm -rf {params.output_dir}/{sra_file}._*.fastq 2>> {log}

                        mv {params.output_dir}/{sra_file}_1.fastq.gz {output.reads[0]} 2>> {log}
                        mv {params.output_dir}/{sra_file}_2.fastq.gz {output.reads[1]} 2>> {log}
                        ''')

                else:
                    r1_list = []
                    r2_list = []
                    for sra in input:
                        sra_file = os.path.basename(sra)
                        r1_list.append(os.path.join(params.output_dir,
                                                 sra_file + "_1.fastq.gz"))
                        r2_list.append(os.path.join(params.output_dir,
                                                 sra_file + "_2.fastq.gz"))
                        shell(
                            f'''
                            rm -rf {params.output_dir}/{sra_file}* 2>> {log}
                            rm -rf {params.output_dir}.{sra_file}.temp 2>> {log}

                            fasterq-dump \
                            --threads {threads} \
                            --split-3 \
                            --temp {params.output_dir}.{sra_file}.temp \
                            --outdir {params.output_dir} {sra} 2>> {log}

                            rm -rf {params.output_dir}.{sra_file}.temp 2>> {log}
                            pigz --processes {threads} {params.output_dir}/{sra_file}_1.fastq 2>> {log}
                            pigz --processes {threads} {params.output_dir}/{sra_file}_2.fastq 2>> {log}
                            rm -rf {params.output_dir}/{sra_file}._*.fastq 2>> {log}
                            ''')

                    r1_str = " ".join(r1_list)
                    r2_str = " ".join(r2_list)
                    shell(f'''cat {r1_str} > {output.reads[0]} 2>> {log}''')
                    shell(f'''cat {r2_str} > {output.reads[1]} 2>> {log}''')
                    shell(f'''rm -rf {r1_str} 2>> {log}''')
                    shell(f'''rm -rf {r2_str} 2>> {log}''')


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
                         "report/stats/{sample}_raw_stats.tsv")
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
                             "report/stats/{sample}_raw_stats.tsv"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "raw_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            metapi.merge(input, metapi.parse, threads, output=output[0])


    rule raw_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "raw_stats.tsv")

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