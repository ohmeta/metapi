rule sra2fq:
    input:
        lambda wildcards: get_sra(_samples, wildcards)
    output:
        r1 = temp(os.path.join(config["results"]["sra2fq"], "{sample}.1.fq.gz")),
        r2 = temp(os.path.join(cofnig["results"]["sra2fq"], "{sample}.2.fq.gz"))
    params:
        sraid = lambda wildcards, input: os.path.basename(input).split(".")[0]
        outdir = config["results"]["sra2fq"],
        reads_direction = str("+"),
        header_format = str("@$ac-$si/$ri")
    shell:
        '''
        fastq-dump --gzip --split-3 \
        --defline-qual '{params.reads_direction}' \
        --defline-seq '{params.header_format}' \
        --outdir {params.outdir} {input}
        mv {params.outdir}/{params.sraid}.1.fastq.gz {output.r1}
        mv {params.outdir}/{params.sraid}.2.fastq.gz {output.r2}
        '''
