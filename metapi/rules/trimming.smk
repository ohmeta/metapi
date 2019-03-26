if config["params"]["trimming"]["oas1"]["do"]:
    rule trimming_oas1:
        input:
            r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
        output:
            r1 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz")),
            r2 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz")),
            single = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.single.fq.gz")),
            stat_out = protected(os.path.join(config["results"]["trimming"], "{sample}.trimmed.stat_out"))
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
            expand(temp("{trimming}/{{sample}}.trimmed.{read}.fq.gz"),
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
            r1 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz")),
            r2 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz")),
            html = protected(os.path.join(config["results"]["trimming"], "{sample}.fastp.html")),
            json = protected(os.path.join(config["results"]["trimming"], "{sample}.fastp.json"))
        params:
            compression = config["params"]["trimming"]["fastp"]["compression"],
            use_slide_window = config["params"]["trimming"]["fastp"]["use_slide_window"],
            cut_front_window_size = config["params"]["trimming"]["fastp"]["cut_front_window_size"],
            cut_front_mean_quality = config["params"]["trimming"]["fastp"]["cut_front_mean_quality"],
            cut_tail_window_size = config["params"]["trimming"]["fastp"]["cut_tail_window_size"],
            cut_tail_mean_quality = config["params"]["trimming"]["fastp"]["cut_tail_mean_quality"],
            cut_right_window_size = config["params"]["trimming"]["fastp"]["cut_right_window_size"],
            cut_right_mean_quality = config["params"]["trimming"]["fastp"]["cut_right_mean_quality"],
            length_required = config["params"]["trimming"]["fastp"]["length_required"],
            n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"],
            adapter_trimming = 'disable_adapter_trimming' if config["params"]["trimming"]["fastp"]["disable_adapter_trimming"] else ""
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.fastp.log")
        threads:
            config["params"]["trimming"]["fastp"]["threads"]
        run:
            if params.use_slide_window:
                shell(
                    '''
                    fastp \
                    --in1 {input.r1} \
                    --in2 {input.r2} \
                    --out1 {output.r1} \
                    --out2 {output.r2} \
                    --compression {params.compression} \
                    --{params.adapter_trimming} \
                    --cut_front \
                    --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_tail_window_size} \
                    --cut_right_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                    ''')
            else:
                shell(
                    '''
                    fastp \
                    --in1 {input.r1} \
                    --in2 {input.r2} \
                    --out1 {output.r1} \
                    --out2 {output.r2} \
                    --compression {params.compression} \
                    --{params.adapter_trimming} \
                    --cut_front \
                    --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --thread {threads} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                    ''')

    rule multiqc_fastp:
        input:
            expand("{trimming}/{sample}.fastp.json",
                   trimming=config["results"]["trimming"],
                   sample=_samples.index)
        output:
            html = protected(os.path.join(config["results"]["trimming"], "fastp_multiqc_report.html")),
            data_dir = directory(os.path.join(config["results"]["trimming"], "fastp_multiqc_report_data"))
        log:
            os.path.join(config["logs"]["trimming"], "multiqc_fastp.log")
        params:
            outdir = config["results"]["trimming"]
        shell:
            '''
            multiqc --outdir {params.outdir} --title fastp --module fastp {input} 2> {log}
            '''
