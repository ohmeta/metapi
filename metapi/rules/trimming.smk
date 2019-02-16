rule trimming:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        r1 = os.path.join(config["results"]["trimming"], "trimmed.{sample}_1.fq.gz"),
        r2 = os.path.join(config["results"]["trimming"], "trimmed.{sample}_2.fq.gz"),
        html = os.path.join(config["results"]["trimming"], "{sample}_fastp.html"),
        json = os.path.join(config["results"]["trimming"], "{sample}_fastp.json")
    params:
        compression = config["params"]["trimming"]["fastp"]["compression"],
        length_required = config["params"]["trimming"]["fastp"]["length_required"],
        n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"]
    log:
        os.path.join(config["logs"]["trimming"], "{sample}.fastp.log")
    threads:
        config["params"]["trimming"]["fastp"]["threads"]
    shell:
        '''
        fastp \
        --in1 {input.r1} \
        --in2 {input.r2} \
        --out1 {output.r1} \
        --out2 {output.r2} \
        --compression {params.compression} \
        --disable_adapter_trimming \
        --cut_front \
        --cut_tail \
        --n_base_limit {params.n_base_limit} \
        --length_required {params.length_required} \
        --html {output.html} \
        --json {output.json} 2> {log}
        '''

