if config["params"]["trimming"]["oas1"]["do"]:
    rule trimming_oas1:
        input:
            r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
        output:
            r1 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz"),
            single = os.path.join(config["results"]["trimming"], "{sample}.trimmed.single.fq.gz"),
            stat_out = os.path.join(config["results"]["trimming"], "{sample}.trimmed.stat_out")
        params:
            prefix = "{sample}",
            qual_system = config["params"]["trimming"]["oas1"]["qual_system"],
            min_length = config["params"]["trimming"]["oas1"]["min_length"],
            seed_oa = config["params"]["trimming"]["oas1"]["seed_oa"],
            fragment_oa = config["params"]["trimming"]["oas1"]["fragment_oa"]
        shell:
            '''
            OAs1 {input.r1},{input.r2} {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}
            '''

if config["params"]["trimming"]["sickle"]["do"]:
    rule trimming_sickle:
        input:
            r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
        output:
            expand("{trimming}/{{sample}}.trimmed.{read}.fq.gz",
                   trimming=config["results"]["trimming"],
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
            r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
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
