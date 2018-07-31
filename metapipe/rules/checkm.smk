rule checkm_lineage_wf:
    input:
        default = os.path.join(config["logs"]["binning"]["metabat2"], "{sample}.metabat2.done"),
        bins_dir = directory(os.path.join(config["results"]["binning"]["bins"], "{sample}.metabat2_out"))
    output:
        checkm_txt = os.path.join(config["results"]["checkm"]["out"],
                                  "{sample}.checkm.txt"),
        checkm_data_dir = directory(os.path.join(config["results"]["checkm"]["data"],
                                       "{sample}"))
    params:
        txt_dir = directory(config["results"]["checkm"]["out"]),
        data_dir = directory(config["results"]["checkm"]["data"]),
        checkm_env = config["params"]["checkm"]["env"]
    threads:
        config["params"]["checkm"]["threads"]
    shell:
        '''
        set +u; source activate {params.checkm_env}; set -u;
        mkdir -p {params.txt_dir}
        mkdir -p {params.data_dir}
        checkm lineage_wf -f {output.checkm_txt} -t {threads} -x fa {input.bins_dir}/ {output.checkm_data_dir}/
        '''
