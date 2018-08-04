rule individual_assembly:
    input:
        reads = expand("{trim}/{{sample}}_{{unit}}.trimmed.{read}.fq.gz",
                       trim=config["results"]["trim"],
                       read=["1", "2"])
    output:
        os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out/{sample}_{unit}.contigs.fa")
    params:
        min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
        megahit_threads = config["params"]["assembly"]["megahit"]["threads"],
        out_dir = os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out"),
        out_prefix = "{sample}_{unit}"
    shell:
        "megahit -1 {input.reads[0]} -2 {input.reads[1]} "
        "-t {params.megahit_threads} "
        "--min-contig-len {params.min_contig} "
        "--out-dir {params.out_dir} "
        "--out-prefix {params.out_prefix} --continue"
