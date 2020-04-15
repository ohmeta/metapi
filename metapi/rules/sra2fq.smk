rule sra2fq:
    input:
        lambda wildcards: metapi.get_reads(SAMPLES, wildcards, "sra")
    output:
        r1 = os.path.join(config["results"]["sra2fq"], "{sample}.1.fq.gz"),
        r2 = os.path.join(config["results"]["sra2fq"], "{sample}.2.fq.gz")
    params:
        outdir = config["results"]["sra2fq"],
        reads_direction = str("+"),
        header_format = str("@$ac-$si/$ri")
    run:
        length = len(input)
        if length == 1:
            sra_id = os.path.basename(input[0]).split(".")[0]
            shell(
                '''
                fastq-dump --gzip --split-3 \
                --defline-qual '{params.reads_direction}' \
                --defline-seq '{params.header_format}' \
                --outdir {params.outdir} {input[0]}
                mv {params.outdir}/%s_1.fastq.gz {output.r1}
                mv {params.outdir}/%s_2.fastq.gz {output.r2}
                ''' % (sra_id, sra_id))
            shell('''rm -rf {params.outdir}/%s.fastq.gz''' % sra_id)
        else:
            r1_list = []
            r2_list = []
            for sra_file in input:
                sra_id = os.path.basename(sra_file).split(".")[0]
                r1_list.append(os.path.join(params.outdir, sra_id + "_1.fastq.gz"))
                r2_list.append(os.path.join(params.outdir, sra_id + "_2.fastq.gz"))
                shell(
                    '''
                    fastq-dump --gzip --split-3 \
                    --defline-qual '{params.reads_direction}' \
                    --defline-seq '{params.header_format}' \
                    --outdir {params.outdir} %s
                    ''' % sra_file)
                shell('''rm -rf {params.outdir}/%s.fastq.gz''' % sra_id)
            r1_str = " ".join(r1_list)
            r2_str = " ".join(r2_list)
            shell('''cat %s > %s''' % (r1_str, output.r1))
            shell('''cat %s > %s''' % (r2_str, output.r2))
            shell('''rm -rf %s''' % r1_str)
            shell('''rm -rf %s''' % r2_str)
