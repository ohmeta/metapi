rule checkm_lineage_wf:
    input:
        default = os.path.join(config["logs"]["binning"]["metabat2"], "{sample}.metabat2.done")
    output:
        checkm_txt = os.path.join(config["results"]["checkm"],
                                  "checkm_out/{sample}.checkm.txt"),
        checkm_data_dir = os.path.join(config["results"]["checkm"],
                                       "checkm_data/{sample}.checkm_out/")
    params:
        bins_dir = config["results"]["binning"]["bins"],
        bins_link_dir = os.path.join(config["results"]["checkm"], "checkm_input"),
        txt_dir = os.path.join(config["results"]["checkm"], "checkm_out"),
        data_dir = os.path.join(config["results"]["checkm"], "checkm_data"),
        checkm_env = config["params"]["checkm"]["env"],
        lineage_threads = config["params"]["checkm"]["threads"]
    shell:
        '''
        set +u; source activate {params.checkm_env}; set -u;
        mkdir -p {params.bins_link_dir}
        mkdir -p {params.txt_dir}
        mkdir -p {params.data_dir}
        find {params.bins_dir} -type f -name "*.bin*" | xargs -I {} ln -s {} {params.bins_link_dir}/
        checkm lineage_wf -f {output.checkm_txt} -t {params.lineage_threads} \
        -x fa {params.bins_link_dir}/ {output.checkm_data_dir}
        '''

'''
rule checkm_filter_wf:
    input:
        checkm_txt = os.path.join(config["results"]["checkm"],
                                  "checkm_out/{sample}.checkm.txt")
    output:
        checkm_filter_txt = os.path.join(config["results"]["checkm"],
                                         "filter/checkm.filter.txt")
    params:
        fa_dir = os.path.join(config["results"]["binning"],
                              "bins/{sample}.metabat2_out"),
        link_dir = os.path.join(config["results"]["checkm"], "checkm_bins/"),
        completeness = config["params"]["checkm"]["completeness"],
        contamination = config["params"]["checkm"]["contamination"]
    run:
        import os
        import re
        import subprocess
        with open(input.checkm_txt, 'r') as checkm_out_h:
            next(checkm_out_h)
            for info in checkm_out_h:
                bin_info = re.split(r'\s+', info.strip())
                if (float(bin_info[-2]) < float(params.contamination)) and (float(bin_info[-3]) > float(params.contamination)):
                    bin_path = os.path.join(params.fa_dir, bin_info[0] + ".fa")
                    cmd = "ln -s %s %s".format(bin_path, params.link_dir)
                    subprocess.run(cmd, shell=true)
'''
