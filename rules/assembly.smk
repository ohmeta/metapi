rule individual_assembly:
    input:
        reads = expand("{rmhost_dir}/{{sample}}_{{unit}}.rmhost.{read}.fq.gz",
                       rmhost_dir=config["results"]["rmhost"],
                       read=["1", "2"])
    output:
        os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out/{sample}_{unit}.contigs.fa")
    params:
        megahit_threads = config["params"]["assembly"]["megahit_threads"],
        out_dir = os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out"),
        out_prefix = "{sample}_{unit}"
    shell:
        "megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {params.megahit_threads} "
        "--out-dir {params.out_dir} --out-prefix {params.out_prefix} "

#rule co_assembly:
#    input:

#    output:
