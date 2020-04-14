def raw_reads(wildcards):
    if IS_PE:
        return [metapi.manager.get_reads(SAMPLES, wildcards, "fq1"),
                metapi.manager.get_reads(SAMPLES, wildcards, "fq2")]
    else:
        return [metapi.manager.get_reads(SAMPLES, wildcards, "fq1")]


if config["params"]["trimming"]["oas1"]["do"]:
    rule trimming_oas1:
        input:
            unpack(raw_reads)
        output:
            reads = temp(expand(os.path.join(config["results"]["trimming"], "{{sample}}.trimmed{read}.fq.gz"),
                                read=[".1", ".2", ".single"] if IS_PE else "")),
            stat_out = os.path.join(config["results"]["trimming"], "{sample}.trimmed.stat_out")
        params:
            prefix = "{sample}",
            qual_system = config["params"]["trimming"]["oas1"]["qual_system"],
            min_length = config["params"]["trimming"]["oas1"]["min_length"],
            seed_oa = config["params"]["trimming"]["oas1"]["seed_oa"],
            fragment_oa = config["params"]["trimming"]["oas1"]["fragment_oa"]
        run:
            reads_num = len(input)
            if IS_PE:
                if reads_num == 2:
                    shell("OAs1 {input[0]},{input[1]} {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}")
                else:
                    r1_str = " ".join(input[0:reads_num//2])
                    r2_str = " ".join(input[reads_num//2:])
                    r1 = os.path.join(config["results"]["trimming"], "%s.raw.1.fq.gz" % params.prefix)
                    r2 = os.path.join(config["results"]["trimming"], "%s.raw.2.fq.gz" % params.prefix)
                    shell("cat %s > %s" % (r1_str, r1))
                    shell("cat %s > %s" % (r2_str, r2))
                    shell("OAs1 %s,%s {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}" % (r1, r2))
                    shell("rm -rf %s" % r1)
                    shell("rm -rf %s" % r2)
            else:
                if reads_num == 1:
                    shell("OAs1 {input[0]} {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}")
                else:
                    r_str = " ".join(input)
                    r = os.path.join(config["results"]["trimming"], "%s.raw.fq.gz" % params.prefix)
                    shell("cat %s > %s" % (r_str, r))
                    shell("OAs1 %s {params.prefix} {params.qual_system} {params.min_length} {params.seed_oa} {params.fragment_oa}" % r)
                    shell("rm -rf %s" % r)


if config["params"]["trimming"]["sickle"]["do"]:
    rule trimming_sickle:
        input:
            unpack(raw_reads)
        output:
            reads = temp(expand(os.path.join(config["results"]["trimming"], "{{sample}}.trimmed{read}.fq.gz"),
                                read=[".1", ".2", ".single"] if IS_PE else ""))
        params:
            qual_type = config["params"]["trimming"]["sickle"]["qual_type"],
            qual_cutoff = config["params"]["trimming"]["sickle"]["qual_cutoff"],
            length_cutoff = config["params"]["trimming"]["sickle"]["length_cutoff"]
        log:
            os.path.join(config["logs"]["trimming"], "{sample}.sickle.log")
        run:
            reads_num = len(input)
            if IS_PE:
                if reads_num == 2:
                    shell("sickle pe --pe-file1 {input[0]} --pe-file2 {input[1]} \
                           --output-pe1 {output.reads[0]} --output-pe2 {output.reads[1]} --output-single {output.reads[2]} \
                           --gzip-output --qual-type {params.qual_type} \
                           --qual-threshold {params.qual_cutoff} --length-threshold {params.length_cutoff} 2> {log}")
                else:
                    r1_str = " ".join(input[0:reads_num//2])
                    r2_str = " ".join(input[reads_num//2:])
                    r1 = os.path.join(config["results"]["trimming"], "%s.raw.1.fq.gz" % params.prefix)
                    r2 = os.path.join(config["results"]["trimming"], "%s.raw.2.fq.gz" % params.prefix)
                    shell("cat %s > %s" % (r1_str, r1))
                    shell("cat %s > %s" % (r2_str, r2))
                    shell("sickle pe --pe-file1 %s --pe-file2 %s \
                           --output-pe1 {output.reads[0]} --output-pe2 {output.reads[1]} --output-single {output.reads[2]} \
                           --gzip-output --qual-type {params.qual_type} \
                           --qual-threshold {params.qual_cutoff} \
                           --length-threshold {params.length_cutoff} 2> {log}" % (r1, r2))
                    shell("rm -rf %s" % r1)
                    shell("rm -rf %s" % r2)
            else:
                if reads_num == 1:
                    shell("sickle se --fastq-file {input[0]} \
                           --output-file {output.reads[0]}\
                           --gzip-output --qual-type {params.qual_type} \
                           --qual-threshold {params.qual_cutoff} --length-threshold {params.length_cutoff} 2> {log}")
                else:
                    r_str = " ".join(input)
                    r = os.path.join(config["results"]["trimming"], "%s.raw.fq.gz" % params.prefix)
                    shell("cat %s > %s" % (r_str, r))
                    shell("sickle se --fastq-file %s \
                           --output-file {output.reads[0]}\
                           --gzip-output --qual-type {params.qual_type} \
                           --qual-threshold {params.qual_cutoff} \
                           --length-threshold {params.length_cutoff} 2> {log}" % r)
                    shell("rm -rf %s" % r)


if config["params"]["trimming"]["fastp"]["do"]:
    rule trimming_fastp:
        input:
            unpack(raw_reads)
        output:
            reads = temp(expand(os.path.join(config["results"]["trimming"], "{{sample}}.trimmed{read}.fq.gz"),
                                read=[".1", ".2"] if IS_PE else "")),
            html = os.path.join(config["results"]["trimming"], "{sample}.fastp.html"),
            json = os.path.join(config["results"]["trimming"], "{sample}.fastp.json")
        params:
            prefix = "{sample}",
            compression = config["params"]["trimming"]["fastp"]["compression"],
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
        run:
            reads_num = len(input)
            if IS_PE:
                if reads_num == 2:
                    if config["params"]["trimming"]["fastp"]["use_slide_window"]:
                        shell("fastp --in1 {input[0]} --in2 {input[1]} --out1 {output.reads[0]} --out2 {output.reads[1]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_right \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_right_window_size {params.cut_right_window_size} \
                               --cut_right_mean_quality {params.cut_right_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}")
                    else:
                        shell("fastp --in1 {input[0]} --in2 {input[1]} --out1 {output.reads[0]} --out2 {output.reads[1]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_tail \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_tail_window_size {params.cut_tail_window_size} \
                               --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}")
                else:
                    r1_str = " ".join(input[0:reads_num//2])
                    r2_str = " ".join(input[reads_num//2:])
                    r1 = os.path.join(config["results"]["trimming"], "%s.raw.1.fq.gz" % params.prefix)
                    r2 = os.path.join(config["results"]["trimming"], "%s.raw.2.fq.gz" % params.prefix)
                    shell("cat %s > %s" % (r1_str, r1))
                    shell("cat %s > %s" % (r2_str, r2))
                    if config["params"]["trimming"]["fastp"]["use_slide_window"]:
                        shell("fastp --in1 %s --in2 %s --out1 {output.reads[0]} --out2 {output.reads[1]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_right \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_right_window_size {params.cut_right_window_size} \
                               --cut_right_mean_quality {params.cut_right_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}" % (r1, r2))
                    else:
                        shell("fastp --in1 %s --in2 %s --out1 {output.reads[0]} --out2 {output.reads[1]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_tail \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_tail_window_size {params.cut_tail_window_size} \
                               --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}" % (r1, r2))
                    shell("rm -rf %s" % r1)
                    shell("rm -rf %s" % r2)
            else:
                if reads_num == 1:
                    if config["params"]["trimming"]["fastp"]["use_slide_window"]:
                        shell("fastp --in1 {input[0]} --out1 {output.reads[0]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_right \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_right_window_size {params.cut_right_window_size} \
                               --cut_right_mean_quality {params.cut_right_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}")
                    else:
                        shell("fastp --in1 {input[0]} --out1 {output.reads[0]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_tail \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_tail_window_size {params.cut_tail_window_size} \
                               --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}")
                else:
                    r_str = " ".join(input)
                    r = os.path.join(config["results"]["trimming"], "%s.raw.fq.gz" % params.prefix)
                    shell("cat %s > %s" % (r_str, r))
                    if config["params"]["trimming"]["fastp"]["use_slide_window"]:
                        shell("fastp --in1 %s --out1 {output.reads[0]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_right \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_right_window_size {params.cut_right_window_size} \
                               --cut_right_mean_quality {params.cut_right_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}" % r)
                    else:
                        shell("fastp --in1 %s --out1 {output.reads[0]} \
                               --compression {params.compression} {params.adapter_trimming} \
                               --cut_front --cut_tail \
                               --cut_front_window_size {params.cut_front_window_size} \
                               --cut_front_mean_quality {params.cut_front_mean_quality} \
                               --cut_tail_window_size {params.cut_tail_window_size} \
                               --cut_tail_mean_quality {params.cut_tail_mean_quality} \
                               --n_base_limit {params.n_base_limit} --length_required {params.length_required} \
                               --thread {threads} --html {output.html} --json {output.json} 2> {log}" % r)
                    shell("rm -rf %s" % r)

    rule multiqc_fastp:
        input:
            expand("{trimming}/{sample}.fastp.json",
                   trimming=config["results"]["trimming"],
                   sample=SAMPLES.index.unique())
        output:
            html = os.path.join(config["results"]["trimming"], "fastp_multiqc_report.html"),
            data_dir = directory(os.path.join(config["results"]["trimming"], "fastp_multiqc_report_data"))
        log:
            os.path.join(config["logs"]["trimming"], "multiqc_fastp.log")
        params:
            outdir = config["results"]["trimming"]
        shell:
            '''
            multiqc --outdir {params.outdir} --title fastp --module fastp {input} 2> {log}
            '''


rule trimming_report:
    input:
        reads = expand(os.path.join(config["results"]["trimming"], "{{sample}}.trimmed{read}.fq.gz"),
                       read=[".1", ".2"] if IS_PE else "")
    output:
        os.path.join(config["results"]["report"]["trimming"], "{sample}.trimming.stats.tsv")
    params:
        fq_encoding = config["params"]["report"]["seqkit"]["fq_encoding"],
        sample_id = "{sample}"
    threads:
        config["params"]["report"]["seqkit"]["threads"]
    run:
        if IS_PE:
            shell("seqkit stats --all --basename --tabular \
                  --fq-encoding %s \
                  --out-file %s \
                  --threads %d %s" % (params.fq_encoding, output, threads, " ".join(input)))
            metapi.qcer.change(output[0], params.sample_id, "trimming", "pe", ["fq1", "fq2"])
        else:
            shell("seqkit stats --all --basename --tabular \
                  --fq-encoding %s \
                  --out-file %s \
                  --threads %d %s" % (params.fq_encoding, output, threads, input))
            metapi.qcer.change(output[0], params.sample_id, "trimming", "se", ["fq1"])


rule merge_trimming_report:
    input:
        expand("{reportout}/{sample}.trimming.stats.tsv",
               reportout=config["results"]["report"]["trimming"],
               sample=SAMPLES.index.unique())
    output:
        os.path.join(config["results"]["report"]["base_dir"], "trimming.stats.tsv")
    run:
        metapi.tooler.merge(input, metapi.qcer.parse, 8, save=True, output=output[0])
