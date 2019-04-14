if config["params"]["trimming"]["oas1"]["do"]:
    rule trimming_oas1:
        input:
            r1 = lambda wildcards: get_reads(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_reads(_samples, wildcards, "fq2")
        output:
            r1 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz")),
            r2 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz")),
            single = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.single.fq.gz")),
            stat_out = protected(os.path.join(config["results"]["trimming"], "{sample}.trimmed.stat_out"))
        params:
            prefix = "{sample}",
            count = lambda wildcards, input: len(input.r1),
            r1_str = lambda wildcards, input: " ".join(input.r1),
            r2_str = lambda wildcards, input: " ".join(input.r2),
            r1 = os.path.join(config["results"]["trimming"], "{sample}.raw.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.raw.2.fq.gz"),
            qual_system = config["params"]["trimming"]["oas1"]["qual_system"],
            min_length = config["params"]["trimming"]["oas1"]["min_length"],
            seed_oa = config["params"]["trimming"]["oas1"]["seed_oa"],
            fragment_oa = config["params"]["trimming"]["oas1"]["fragment_oa"]
        shell:
            '''
            if [ {params.count} -eq 1 ]; then
                OAs1 {input.r1[0]},{input.r2[0]} {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}
            else
                cat {params.r1_str} > {params.r1}
                cat {params.r2_str} > {params.r2}
                OAs1 {params.r1},{params.r2} {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}
                rm -rf {params.r1}
                rm -rf {params.r2}
            fi
            '''


if config["params"]["trimming"]["sickle"]["do"]:
    rule trimming_sickle:
        input:
            r1 = lambda wildcards: get_reads(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_reads(_samples, wildcards, "fq2")
        output:
            expand(temp("{trimming}/{{sample}}.trimmed.{read}.fq.gz"),
                   trimming=config["results"]["trimming"],
                   read=["1", "2", "single"])
        params:
            count = lambda wildcards, input: len(input.r1),
            r1_str = lambda wildcards, input: " ".join(input.r1),
            r2_str = lambda wildcards, input: " ".join(input.r2),
            r1 = os.path.join(config["results"]["trimming"], "{sample}.raw.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.raw.2.fq.gz"),
            qual_type = config["params"]["trimming"]["sickle"]["qual_type"],
            qual_cutoff = config["params"]["trimming"]["sickle"]["qual_cutoff"],
            length_cutoff = config["params"]["trimming"]["sickle"]["length_cutoff"]
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.sickle.log")
        shell:
            '''
            if [ {params.count} -eq 1 ]; then
                sickle pe -f {input.r1[0]} -r {input.r2[0]} \
                -o {output[0]} -p {output[1]} -s {output[2]} \
                --gzip-output \
                -t {params.qual_type} \
                -q {params.qual_cutoff} \
                -l {params.length_cutoff} 2> {log}
            else
                cat {params.r1_str} > {params.r1}
                cat {params.r2_str} > {params.r2}
                sickle pe -f {params.r1} -r {params.r2} \
                -o {output[0]} -p {output[1]} -s {output[2]} \
                --gzip-output \
                -t {params.qual_type} \
                -q {params.qual_cutoff} \
                -l {params.length_cutoff} 2> {log}
                rm -rf {params.r1}
                rm -rf {params.r2}
            fi
            '''


if config["params"]["trimming"]["fastp"]["do"]:
    rule trimming_fastp:
        input:
            r1 = lambda wildcards: get_reads(_samples, wildcards, "fq1"),
            r2 = lambda wildcards: get_reads(_samples, wildcards, "fq2")
        output:
            r1 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.1.fq.gz")),
            r2 = temp(os.path.join(config["results"]["trimming"], "{sample}.trimmed.2.fq.gz")),
            html = protected(os.path.join(config["results"]["trimming"], "{sample}.fastp.html")),
            json = protected(os.path.join(config["results"]["trimming"], "{sample}.fastp.json"))
        params:
            count = lambda wildcards, input: len(input.r1),
            r1_str = lambda wildcards, input: " ".join(input.r1),
            r2_str = lambda wildcards, input: " ".join(input.r2),
            r1 = os.path.join(config["results"]["trimming"], "{sample}.raw.1.fq.gz"),
            r2 = os.path.join(config["results"]["trimming"], "{sample}.raw.2.fq.gz"),
            compression = config["params"]["trimming"]["fastp"]["compression"],
            use_slide_window = "true" if config["params"]["trimming"]["fastp"]["use_slide_window"] else "false",
            cut_front_window_size = config["params"]["trimming"]["fastp"]["cut_front_window_size"],
            cut_front_mean_quality = config["params"]["trimming"]["fastp"]["cut_front_mean_quality"],
            cut_tail_window_size = config["params"]["trimming"]["fastp"]["cut_tail_window_size"],
            cut_tail_mean_quality = config["params"]["trimming"]["fastp"]["cut_tail_mean_quality"],
            cut_right_window_size = config["params"]["trimming"]["fastp"]["cut_right_window_size"],
            cut_right_mean_quality = config["params"]["trimming"]["fastp"]["cut_right_mean_quality"],
            length_required = config["params"]["trimming"]["fastp"]["length_required"],
            n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"],
            adapter_trimming = '--disable_adapter_trimming' if config["params"]["trimming"]["fastp"]["disable_adapter_trimming"] else ""
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.fastp.log")
        threads:
            config["params"]["trimming"]["fastp"]["threads"]
        shell:
            '''
            if [ {params.count} -gt 1 ]; then
                cat {params.r1_str} > {params.r1}
                cat {params.r2_str} > {params.r2}
                if {params.use_slide_window}; then
                    fastp --in1 {params.r1} --in2 {params.r2} --out1 {output.r1} --out2 {output.r2} \
                    --compression {params.compression} {params.adapter_trimming} \
                    --cut_front --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_tail_window_size} \
                    --cut_right_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                    --thread {threads} --html {output.html} --json {output.json} 2> {log}
                else
                    fastp --in1 {params.r1} --in2 {params.r2} --out1 {output.r1} --out2 {output.r2} \
                    --compression {params.compression} {params.adapter_trimming} \
                    --cut_front --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                    --thread {threads} --html {output.html} --json {output.json} 2> {log}
                fi
                rm -rf {params.r1}
                rm -rf {params.r2}
            else
                if {params.use_slide_window}; then
                    fastp --in1 {input.r1[0]} --in2 {input.r2[0]} --out1 {output.r1} --out2 {output.r2} \
                    --compression {params.compression} {params.adapter_trimming} \
                    --cut_front --cut_right \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_right_window_size {params.cut_tail_window_size} \
                    --cut_right_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                    --thread {threads} --html {output.html} --json {output.json} 2> {log}
                else
                    fastp --in1 {input.r1[0]} --in2 {input.r2[0]} --out1 {output.r1} --out2 {output.r2} \
                    --compression {params.compression} {params.adapter_trimming} \
                    --cut_front --cut_tail \
                    --cut_front_window_size {params.cut_front_window_size} \
                    --cut_front_mean_quality {params.cut_front_mean_quality} \
                    --cut_tail_window_size {params.cut_tail_window_size} \
                    --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                    --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                    --thread {threads} --html {output.html} --json {output.json} 2> {log}
                fi
            fi
            '''


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
