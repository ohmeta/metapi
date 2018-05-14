rule one_assembly:
    input:
        reads = expand("{rmhost_dir}/{{sample}}_{{unit}}.rmhost.{read}.fq.gz",
                       rmhost_dir=config["results"]["rmhost"],
                       read=["1", "2"])
    output:
        os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out/{sample}_{unit}.contigs.fa")
    threads:
        config["params"]["assembly"]["megahit_threads"]
    params:
        out_dir = os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out"),
        out_prefix = "{sample}_{unit}"
    shell:
        "megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} "
        "--out-dir {params.out_dir} --out-prefix {params.out_prefix} --continue"

#rule co_assembly:
#    input:

#    output:
