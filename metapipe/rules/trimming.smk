if config["params"]["trimming"]["sickle"]["do"]:
    rule trimming_sickle:
        input:
            r1 = lambda wildcards: _get_raw_fastq(wildcards, "fq1"),
            r2 = lambda wildcards: _get_raw_fastq(wildcards, "fq2")
        output:
            expand("{trim}/{{sample}}.trimmed.{read}.fq.gz",
                   trim=config["results"]["trimming"],
                   read=["1", "2", "single"])
        params:
            qual_type = config["params"]["trimming"]["sickle"]["qual_type"],
            qual_cutoff = config["params"]["trimming"]["sickle"]["qual_cutoff"],
            length_cutoff = config["params"]["trimming"]["sickle"]["length_cutoff"]
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.sickle.log")
        shell:
            '''
            sickle pe -f {input.r1} -r {input.r2} \
            -o {output[0]} -p {output[1]} -s {output[2]} \
            --gzip-output \
            -t {params.qual_type} \
            -q {params.qual_cutoff} \
            -l {params.length_cutoff} 2> {log}
            '''

if config["params"]["trimming"]["fastp"]["do"]:
    rule trimming_fastp:
        input:
            r1 = lambda wildcards: _get_raw_fastq(wildcards, "fq1"),
            r2 = lambda wildcards: _get_raw_fastq(wildcards, "fq2")
        output:
            r1 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz"),
            html = os.path.join(config["results"]["trimming"], "{sample}.fastp.html"),
            json = os.path.join(config["results"]["trimming"], "{sample}.fastp.json")
        params:
            compression = config["params"]["trimming"]["fastp"]["compression"],
            cut_mean_quality = config["params"]["trimming"]["fastp"]["cut_mean_quality"],
            length_required = config["params"]["trimming"]["fastp"]["length_required"]
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.fastp.log")
        threads:
            config["params"]["trimming"]["fastp"]["threads"]
        shell:
            '''
            fastp --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --compression {params.compression} \
            --cut_mean_quality {params.cut_mean_quality} \
            --length_required {params.length_required} \
            --html {output.html} \
            --json {output.json} 2> {log}
            '''
