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
            use_slide_window = config["params"]["trimming"]["fastp"]["use_slide_window"],
            cut_front_window_size = config["params"]["trimming"]["fastp"]["cut_front_window_size"],
            cut_front_mean_quality = config["params"]["trimming"]["fastp"]["cut_front_mean_quality"],
            cut_tail_window_size = config["params"]["trimming"]["fastp"]["cut_tail_window_size"],
            cut_tail_mean_quality = config["params"]["trimming"]["fastp"]["cut_tail_mean_quality"],
            cut_right_window_size = config["params"]["trimming"]["fastp"]["cut_right_window_size"],
            cut_right_mean_quality = config["params"]["trimming"]["fastp"]["cut_right_mean_quality"],
            length_required = config["params"]["trimming"]["fastp"]["length_required"],
            n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"]
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.fastp.log")
        threads:
            config["params"]["trimming"]["fastp"]["threads"]
        run:
            if {params.use_slide_window}:
                print({input.r1})
                print({input.r2})
                shell(
                    '''
                    fastp \
                    --in1 {input.r1} \
                    --in2 {input.r2} \
                    --out1 {output.r1} \
                    --out2 {output.r2} \
                    --compression {params.compression} \
                    --disable_adapter_trimming \
                    --cut_front \
                    --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_tail_window_size} \
                    --cut_right_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
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
                    --disable_adapter_trimming \
                    --cut_front \
                    --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} \
                    --length_required {params.length_required} \
                    --html {output.html} \
                    --json {output.json} 2> {log}
                    ''')

