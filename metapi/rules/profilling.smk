rule metaphlan2_profilling:
    input:
        reads = clean_reads
    output:
        bt2_out = os.path.join(config["results"]["profilling"]["metaphlan2"]["bowtie2_out"], "{sample}.bowtie2.gz"),
        profile = os.path.join(config["results"]["profilling"]["metaphlan2"]["profile"], "{sample}.metaphlan2.profile")
    params:
        metaphlan2_env = config["params"]["profilling"]["metaphlan2"]["env"],
        input_type = config["params"]["profilling"]["metaphlan2"]["input_type"],
    log:
        os.path.join(config["logs"]["profilling"]["metaphlan2"], "{sample}_metaphlan2.log")
    threads:
        config["params"]["profilling"]["metaphlan2"]["threads"]
    shell:
        '''
        set +u; source activate {params.metaphlan2_env}; set -u;
        metaphlan2.py {input.reads[0]},{input.reads[1]} --bowtie2out {output.bt2_out} --nproc {threads} --input_type {params.input_type} > {output.profile} 2> {log}
        '''

rule metaphlan2_merge:
    input:
        expand("{profile}/{sample}.metaphlan2.profile",
               profile=config["results"]["profilling"]["metaphlan2"]["profile"],
               sample=_samples.index)
    output:
        os.path.join(config["results"]["profilling"]["metaphlan2"]["base_dir"], "metaphlan2.merged.profile")
    log:
        os.path.join(config["logs"]["profilling"]["metaphlan2"], "metaphlan2.merged.log")
    params:
        metaphlan2_env = config["params"]["profilling"]["metaphlan2"]["env"]
    shell:
        '''
        set +u; source activate {params.metaphlan2_env}; set -u;
        merge_metaphlan_tables.py {input} > {output} 2> {log}
        '''
