rule metaphlan2_profilling:
    input:
        reads = clean_reads
    output:
        bt2_out = os.path.join(config["results"]["profilling"]["metaphlan2"]["bowtie2_out"], "{sample}.bowtie2.bz2"),
        profile = os.path.join(config["results"]["profilling"]["metaphlan2"]["profile"], "{sample}.metaphlan2.profile")
    params:
        metaphlan2_env = config["params"]["profilling"]["metaphlan2"]["env"],
        metaphlan2_script = config["params"]["profilling"]["metaphlan2"]["script"],
        bowtie2_exe = config["params"]["profilling"]["metaphlan2"]["bowtie2_exe"],
        input_type = config["params"]["profilling"]["metaphlan2"]["input_type"],
        mpa_pkl = config["params"]["profilling"]["metaphlan2"]["mpa_pkl"],
        bowtie2db = config["params"]["profilling"]["metaphlan2"]["bowtie2db"],
        min_cu_len = config["params"]["profilling"]["metaphlan2"]["min_cu_len"],
        sample_id = "{sample}"
    log:
        os.path.join(config["logs"]["profilling"]["metaphlan2"], "{sample}_metaphlan2.log")
    threads:
        config["params"]["profilling"]["metaphlan2"]["threads"]
    shell:
        '''
        #set +u; source activate {params.metaphlan2_env}; set -u;
        /ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv2/bin/python \
        {params.metaphlan2_script} \
        {input.reads[0]},{input.reads[1]} \
        --bowtie2_exe {params.bowtie2_exe} \
        --mpa_pkl {params.mpa_pkl} \
        --bowtie2db {params.bowtie2db} \
        --min_cu_len {params.min_cu_len} \
        --bowtie2out {output.bt2_out} \
        --nproc {threads} \
        --input_type {params.input_type} \
        --sample_id {params.sample_id} \
        > {output.profile} \
        2> {log}
        '''

rule metaphlan2_merge:
    input:
        expand("{profile}/{sample}.metaphlan2.profile",
               profile=config["results"]["profilling"]["metaphlan2"]["profile"],
               sample=_samples.index.unique())
    output:
        os.path.join(config["results"]["profilling"]["metaphlan2"]["base_dir"], "metaphlan2.merged.profile")
    log:
        os.path.join(config["logs"]["profilling"]["metaphlan2"], "metaphlan2.merged.log")
    params:
        metaphlan2_env = config["params"]["profilling"]["metaphlan2"]["env"]
    shell:
        '''
        #!/bin/bash
        #set +u; source activate {params.metaphlan2_env}; set -u;
        /ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv2/bin/python /ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv2/bin/merge_metaphlan_tables.py {input} > {output} 2> {log}
        sed -i 's/.metaphlan2//g' {output}
        '''
