import os
from find_path import find_path_tag

(R1, R2, RS, RT) = find_path_tag(config["filter_reads_dir"], "clean")

rule raw_report:
    input:
        r1 =
        r2 =
        stat_out =
    output:
    shell:

rule filter_report:
    input:
        r1       = lambda wildcards: R1[wildcards.sample]
        r2       = lambda wildcards: R2[wildcards.sample]
        single   = lambda wildcards: RS[wildcards.sample]
        stat_out = lambda wildcards: RT[wildcards.sample]
    output:
        expand("{filter_report_dir}/01.filter/{sample}.csv",
               filter_report_dir=config["qc_report_dir"])
    params:
        filter_reporter=config["filter_report"]
        stat_out_dir=config["filter_reads_dir"]
        out_dir=config["qc_report_dir"]
    shell:
        "python {params.filter_reporter} -d {params.stat_out_dir} -n 30 -m 100"

rule rmhost_report:
    input:
        r1 =
        r2 =
        single =
        stat_out =
    output:
        expand("{rmhost_report_dir}/02.rmhost/{sample}.tsv",
               rmhost_report_dir=config["qc_report_dir"])
    shell:
        "python {params.rmhost_reporter} -1 input.r1"

rule merger:
    input:
        expand("{report_dir}/clean/{sample}.csv", report_dir=config["qc_report_dir"])
    output:

