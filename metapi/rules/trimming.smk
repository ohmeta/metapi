def raw_samples(wildcards):
    if config["params"]["simulation"]["do"]:
        return [metapi.get_reads(SAMPLES, wildcards, "fq1")[0],
                metapi.get_reads(SAMPLES, wildcards, "fq2")[0]]
    else:
        if IS_PE:
            return [metapi.get_reads(SAMPLES, wildcards, "fq1"),
                    metapi.get_reads(SAMPLES, wildcards, "fq2")]
        else:
            return [metapi.get_reads(SAMPLES, wildcards, "fq1")]


def raw_reads(wildcards, have_single):
    if have_single:
        return expand(os.path.join("{short_reads}",
                                   "{sample}.link_or_merge.out/{sample}.raw{read}.fq.gz"),
                      short_reads=config["results"]["raw"]["short_reads"],
                      read=[".1", ".2", ".single"] if IS_PE else "",
                      sample=wildcards.sample)
    else:
        return expand(os.path.join("{short_reads}",
                                   "{sample}.link_or_merge.out/{sample}.raw{read}.fq.gz"),
                      short_reads=config["results"]["raw"]["short_reads"],
                      read=[".1", ".2"] if IS_PE else "",
                      sample=wildcards.sample)


rule prepare_reads:
    input:
        unpack(raw_samples)
    output:
        expand(os.path.join("{short_reads}",
                            "{{sample}}.link_or_merge.out/{{sample}}.raw{read}.fq.gz"),
               short_reads=config["results"]["raw"]["short_reads"],
               read=[".1", ".2"] if IS_PE else "")
    params:
        output_dir = os.path.join(config["results"]["raw"]["short_reads"],
                                  "{sample}.link_or_merge.out"),
        output_prefix = os.path.join(config["results"]["raw"]["short_reads"],
                                  "{sample}.link_or_merge.out/{sample}")
    run:
        shell("mkdir -p {params.output_dir}")
        reads_num = len(input)
        if IS_PE:
            if reads_num == 2:
                shell("ln -s {input[0]}  {output[0]}")
                shell("ln -s {input[1]}  {output[1]}")
            else:
                r1_str = " ".join(input[0:reads_num//2])
                r2_str = " ".join(input[reads_num//2:])
                r1 = "%s.raw.1.fq.gz" % params.output_prefix
                r2 = "%s.raw.2.fq.gz" % params.output_prefix
                shell("cat %s > %s" % (r1_str, r1))
                shell("cat %s > %s" % (r2_str, r2))
        else:
            if reads_num == 1:
                shell("ln -s {input[0]}  {output[0]}")
            else:
                r_str = " ".join(input)
                r = "%s.raw.fq.gz" % params.output_prefix
                shell("cat %s > %s" % (r_str, r))


rule trimming_oas1:
    input:
        lambda wildcards: raw_reads(wildcards, False)
    output:
        reads = expand(os.path.join("{short_reads}",
                                    "{{sample}}.oas1.out/{{sample}}.trimmed{read}.fq.gz"),
                       short_reads=config["results"]["trimmed"]["short_reads"],
                       read=[".1", ".2", ".single"] if IS_PE else ""),
        stat_out = expand(os.path.join("{short_reads}",
                                       "{{sample}}.oas1.out/{{sample}}.trimmed.stat_out"),
                          short_reads=config["results"]["trimmed"]["short_reads"])
    params:
        output_prefix = os.path.join(config["results"]["trimmed"]["short_reads"],
                                     "{sample}.oas1.out/{sample}"),
        short_reads = os.path.join(config["results"]["trimmed"]["short_reads"],
                                   "{sample}.oas1.out"),
        quality_system = config["params"]["trimming"]["oas1"]["quality_system"],
        min_length = config["params"]["trimming"]["oas1"]["min_length"],
        seed_oa = config["params"]["trimming"]["oas1"]["seed_oa"],
        fragment_oa = config["params"]["trimming"]["oas1"]["fragment_oa"]
    log:
        os.path.join(config["logs"]["trimming"], "trimming.oas1.{sample}.log")
    run:
        reads_str = ",".join(input)
        shell("OAs1 %s \
              {params.output_prefix} \
              {params.quality_system} \
              {params.min_length} \
              {params.seed_oa} \
              {params.fragment_oa} 2> {log}" % reads_str)
        if IS_PE:
            shell("mv {params.output_prefix}.clean.1.fq.gz {output.reads[0]}")
            shell("mv {params.output_prefix}.clean.2.fq.gz {output.reads[1]}")
            shell("mv {params.output_prefix}.clean.single.fq.gz {output.reads[2]}")
        else:
            shell("mv {params.output_prefix}.clean.fq.gz  {output.reads[0]}")
        shell("mv {params.output_prefix}.clean.stat_out {output.stat_out}")


rule trimming_sickle:
    input:
        lambda wildcards: raw_reads(wildcards, False)
    output:
        expand(os.path.join("{short_reads}",
                            "{{sample}}.sickle.out/{{sample}}.trimmed{read}.fq.gz"),
               short_reads=config["results"]["trimmed"]["short_reads"],
               read=[".1", ".2", ".single"] if IS_PE else "")
    params:
        output_prefix = os.path.join(config["results"]["trimmed"]["short_reads"],
                                     "{sample}.sickle.out/{sample}"),
        quality_type = config["params"]["trimming"]["sickle"]["quality_type"],
        quality_cutoff = config["params"]["trimming"]["sickle"]["quality_cutoff"],
        length_cutoff = config["params"]["trimming"]["sickle"]["length_cutoff"]
    log:
        os.path.join(config["logs"]["trimming"], "trimming.sickle.{sample}.log")
    run:
        if IS_PE:
            shell("sickle pe \
                  --pe-file1 {input[0]} \
                  --pe-file2 {input[1]} \
                  --output-pe1 {output[0]} \
                  --output-pe2 {output[1]} \
                  --output-single {output[2]} \
                  --qual-type {params.quality_type} \
                  --qual-threshold {params.quality_cutoff} \
                  --length-threshold {params.length_cutoff} \
                  --gzip-output 2> {log}")
        else:
            shell("sickle se \
                  --fastq-file {input[0]} \
                  --output-file {output[0]} \
                  --qual-type {params.quality_type} \
                  --qual-threshold {params.quality_cutoff} \
                  --length-threshold {params.length_cutoff} \
                  --gzip-output 2> {log}")


rule trimming_fastp:
    input:
        reads = lambda wildcards: raw_reads(wildcards, False)
    output:
        reads = expand(os.path.join("{short_reads}",
                                    "{{sample}}.fastp.out/{{sample}}.trimmed{read}.fq.gz"),
                       short_reads=config["results"]["trimmed"]["short_reads"],
                       read=[".1", ".2", ".single"] if IS_PE else ""),
        html = os.path.join(config["results"]["trimmed"]["short_reads"],
                            "{sample}.fastp.out/{sample}.fastp.html"),
        json = os.path.join(config["results"]["trimmed"]["short_reads"],
                            "{sample}.fastp.out/{sample}.fastp.json")
    params:
        output_prefix = os.path.join(config["results"]["trimmed"]["short_reads"],
                                     "{sample}.fastp.out/{sample}"),
        compression = config["params"]["trimming"]["fastp"]["compression"],
        cut_front_window_size = config["params"]["trimming"]["fastp"]["cut_front_window_size"],
        cut_front_mean_quality = config["params"]["trimming"]["fastp"]["cut_front_mean_quality"],
        cut_tail_window_size = config["params"]["trimming"]["fastp"]["cut_tail_window_size"],
        cut_tail_mean_quality = config["params"]["trimming"]["fastp"]["cut_tail_mean_quality"],
        cut_right_window_size = config["params"]["trimming"]["fastp"]["cut_right_window_size"],
        cut_right_mean_quality = config["params"]["trimming"]["fastp"]["cut_right_mean_quality"],
        length_required = config["params"]["trimming"]["fastp"]["length_required"],
        n_base_limit = config["params"]["trimming"]["fastp"]["n_base_limit"],
        adapter_trimming = '--disable_adapter_trimming' \
            if config["params"]["trimming"]["fastp"]["disable_adapter_trimming"] else ""
    log:
        os.path.join(config["logs"]["trimming"], "trimming.fastp.{sample}.log")
    threads:
        config["params"]["trimming"]["fastp"]["threads"]
    run:
        if IS_PE:
            if config["params"]["trimming"]["fastp"]["use_slide_window"]:
                shell("fastp \
                      --in1 {input[0]} \
                      --in2 {input[1]} \
                      --out1 {output.reads[0]} \
                      --out2 {output.reads[1]} \
                      --compression {params.compression} \
                      {params.adapter_trimming} \
                      --cut_front \
                      --cut_right \
                      --cut_front_window_size {params.cut_front_window_size} \
                      --cut_front_mean_quality {params.cut_front_mean_quality} \
                      --cut_right_window_size {params.cut_right_window_size} \
                      --cut_right_mean_quality {params.cut_right_mean_quality} \
                      --n_base_limit {params.n_base_limit} \
                      --length_required {params.length_required} \
                      --thread {threads} \
                      --html {output.html} \
                      --json {output.json} 2> {log}")
            else:
                shell("fastp \
                      --in1 {input[0]} \
                      --in2 {input[1]} \
                      --out1 {output.reads[0]} \
                      --out2 {output.reads[1]} \
                      --compression {params.compression} \
                      {params.adapter_trimming} \
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
                      --json {output.json} 2> {log}")
        else:
            if config["params"]["trimming"]["fastp"]["use_slide_window"]:
                shell("fastp
                      --in1 {input[0]} \
                      --out1 {output.reads[0]} \
                      --compression {params.compression} \
                      {params.adapter_trimming} \
                      --cut_front \
                      --cut_right \
                      --cut_front_window_size {params.cut_front_window_size} \
                      --cut_front_mean_quality {params.cut_front_mean_quality} \
                      --cut_right_window_size {params.cut_right_window_size} \
                      --cut_right_mean_quality {params.cut_right_mean_quality} \
                      --n_base_limit {params.n_base_limit} \
                      --length_required {params.length_required} \
                      --thread {threads} \
                      --html {output.html} \
                      --json {output.json} 2> {log}")
            else:
                shell("fastp \
                      --in1 {input[0]} \
                      --out1 {output.reads[0]} \
                      --compression {params.compression} \
                      {params.adapter_trimming} \
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
                      --json {output.json} 2> {log}")


rule multiqc_fastp:
    input:
        expand(os.path.join("{short_reads}",
                            "{sample}.fastp.out/{sample}.fastp.json"),
               short_reads=config["results"]["trimming"]["short_reads"],
               sample=SAMPLES.index.unique())
    output:
        html = os.path.join(config["results"]["trimming"]["report"],
                            "report_multiqc_fastp.html"),
        data_dir = directory(os.path.join(config["results"]["trimming"]["report"],
                                          "report_multiqc_fastp_data"))
    log:
        os.path.join(config["logs"]["trimming"], "trimming.multiqc.fastp.log")
    params:
        outdir = config["results"]["trimming"]["report"]
    shell:
        '''
        multiqc --outdir {params.outdir} --title fastp --module fastp {input} 2> {log}
        '''


rule trimming_report:
    input:
        reads = expand(os.path.join("{short_reads}",
                                    "{{sample}}.{{trimmer}}.out/{{sample}}.trimmed{read}.fq.gz"),
                       read=[".1", ".2"] if IS_PE else "")
    output:
        os.path.join(config["results"]["trimming"]["report"], "temp/{sample}.{trimmer}.stats.tsv")
    params:
        fq_encoding = config["params"]["qc_report"]"seqkit"]["fq_encoding"],
        sample_id = "{sample}"
    threads:
        config["params"]["qc_report"]["seqkit"]["threads"]
    run:
        if IS_PE:
            shell("seqkit stats \
                  --all \
                  --basename \
                  --tabular \
                  --fq-encoding %s \
                  --out-file %s \
                  --threads %d %s" % (params.fq_encoding, output, threads, " ".join(input)))
            metapi.change(output[0], params.sample_id, "trimming", "pe", ["fq1", "fq2"])
        else:
            shell("seqkit stats \
                  --all \
                  --basename \
                  --tabular \
                  --fq-encoding %s \
                  --out-file %s \
                  --threads %d %s" % (params.fq_encoding, output, threads, input))
            metapi.change(output[0], params.sample_id, "trimming", "se", ["fq1"])


rule trimming_report_merge:
    input:
        expand(os.path.join("{report}",
                            "temp/{sample}.{{trimmer}}.stats.tsv"),
               report=config["results"]["trimming"]["report"],
               sample=SAMPLES.index.unique())
    output:
        os.path.join(config["results"]["trimming"]["report"],
                     "trimming.{trimmer}.stats.tsv")
    run:
        metapi.merge(input, metapi.parse, 8, save=True, output=output[0])


rule all_output_oas1:
    input:
        expand([
            os.path.join("{short_reads}",
                         "{sample}.oas1.out/{sample}.trimmed{read}.fq.gz"),
            os.path.join("{short_reads}",
                         "{sample}.oas1.out/{sample}.trimmed.stat_out")],
            short_reads=config["results"]["trimming"]["short_reads"],
            read=[".1", ".2", ".single"] if IS_PE else "",
            sample=SAMPLES.index.unique())


rule all_output_sickle:
    input:
        expand(
            os.path.join("{short_reads}",
                         "{sample}.sickle.out/{sample}.trimmed{read}.fq.gz"),
            short_reads=config["results"]["trimming"]["short_reads"],
            read=[".1", ".2", ".single"] if IS_PE else "",
            sample=SAMPLES.index.unique())


rule all_output_fastp:
    input:
        expand([
            os.path.join("{short_reads}", "{sample}.fastp.out/{sample}.trimmed{read}.fq.gz"),
            os.path.join("{short_reads}", "{sample}.fastp.out/{sample}.fastp.html"),
            os.path.join("{short_reads}", "{sample}.fastp.out/{sample}.fastp.json"),
            os.path.join("{report}", "report_multiqc_fastp.html"),
            os.path.join("{report}", "report_multiqc_fastp_data")],
            short_reads=config["results"]["trimming"]["short_reads"],
            report=config["results"]["trimming"]["report"],
            read=[".1", ".2"] if IS_PE else "",
            sample=SAMPLES.index.unique())


rule all_output_trimming_report:
    input:
        os.path.join(config["results"]["trimming"]["report"],
                     "trimming.{trimmer}.stats.tsv")


rule all_output_trimming:
    input:
        rules.all_output_oas1.input,
        rules.all_output_sickle.input,
        rules.all_output_fastp.input,
        rules.all_output_trimming_report
