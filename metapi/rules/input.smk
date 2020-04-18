def raw_reads(wildcards):
    if READS_FORMAT == "fastq":
        if config["params"]["simulate"]["do"]:
            return [[metapi.get_reads(SAMPLES, wildcards, "fq1")[0]],
                    [metapi.get_reads(SAMPLES, wildcards, "fq2")[0]]]
        else:
            if IS_PE:
                return [metapi.get_reads(SAMPLES, wildcards, "fq1"),
                        metapi.get_reads(SAMPLES, wildcards, "fq2")]
            else:
                return [metapi.get_reads(SAMPLES, wildcards, "fq1")]
    elif READS_FORMAT == "sra":
        return [metapi.get_reads(SAMPLES, wildcards, "sra")]


rule prepare_reads:
    input:
        unpack(raw_reads)
    output:
        expand(
            os.path.join(config["output"]["raw"],
                         "short_reads/{{sample}}/{{sample}}{read}.fq.gz"),
            read=[".1", ".2"] if IS_PE else "")
    params:
        output_dir = os.path.join(config["output"]["raw"],
                                  "short_reads/{sample}"),
        output_prefix = os.path.join(config["output"]["raw"],
                                     "short_reads/{sample}/{sample}")
    run:
        reads_num = len(input)

        if READS_FORMAT == "fastq":
            if IS_PE:
                if reads_num == 2:
                    os.symlink(os.path.realpath(input[0]), output[0])
                    os.symlink(os.path.realpath(input[1]), output[1])
                else:
                    r1_str = " ".join(input[0:reads_num//2])
                    r2_str = " ".join(input[reads_num//2:])
                    r1 = "%s.raw.1.fq.gz" % params.output_prefix
                    r2 = "%s.raw.2.fq.gz" % params.output_prefix
                    shell('''cat %s > %s''' % (r1_str, r1))
                    shell('''cat %s > %s''' % (r2_str, r2))
            else:
                if reads_num == 1:
                    os.symlink(os.path.realpath(input[0]), output[0])
                else:
                    r_str = " ".join(input)
                    r = "%s.raw.fq.gz" % params.output_prefix
                    shell('''cat %s > %s''' % (r_str, r))

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
                          "short_reads/{sample}/{sample}{read}.fq.gz"),
             read=[".1", ".2"] if IS_PE else "",
             sample=SAMPLES.index.unique())


def get_reads(wildcards, step, have_single):
    if have_single:
        return expand(
            os.path.join(config["output"][step],
                         "short_reads/{sample}/{sample}{read}.fq.gz"),
            read=[".1", ".2", ".single"] if IS_PE else "",
            sample=wildcards.sample)
    else:
        return expand(
            os.path.join(config["output"][step],
                         "short_reads/{sample}/{sample}{read}.fq.gz"),
            read=[".1", ".2"] if IS_PE else "",
            sample=wildcards.sample)
