rule checkm_lineage_wf:
    input:
        fa_dir = os.path.join(config["results"]["binning"], "bins/{sample}_{unit}.metabat2_out/"),
        default = os.path.join(config["logs"]["binning"], "{sample}_{unit}.done")
    output:
        checkm_txt = os.path.join(config["results"]["checkm"], "checkm_out/{sample}_{unit}.checkm.txt"),
        checkm_data_dir = os.path.join(config["results"]["checkm"], "checkm_data/{sample}_{unit}.checkm_out")
    params:
        data_dir = os.path.join(config["results"]["checkm"], "checkm_data"),
        checkm_env = config["params"]["checkm"]["env"],
        lineage_threads = config["params"]["checkm"]["threads"]
    shell:
        '''
        set +u; source activate {params.checkm_env}; set -u;
        rm -rf {output.checkm_data_dir}
        checkm lineage_wf -f {output.checkm_txt} -t {params.lineage_threads} \
        -x fa {input.fa_dir} {output.checkm_data_dir} -r
        '''
'''
rule checkm_filter_wf:
    input:
        checkm_txt = os.path.join(config["results"]["checkm"], "checkm_out/{sample}_{unit}.checkm.txt")
    output:
        filter_txt = os.path.join(config["results"]["checkm"], "checkm.filter.txt")
    params:
        fa_dir = os.path.join(config["results"]["binning"], "bins/{sample}_{unit}.metabat2_out"),
        link_dir = os.path.join(config["results"]["checkm"],"checkm_bins/"),
        completeness = config["params"]["checkm"]["completeness"],
        contamination = config["params"]["checkm"]["contamination"]
    run:
        import os
        import re
        import subprocess
        with open(input.checkm_txt, 'r') as checkm_h, open(output.filter_txt, 'aw') as filter_h:
            next(checkm_h)
            for info in checkm_h:
                bin_info = re.split(r'\s+', info.strip())
                if (float(bin_info[-2]) < float(params.contamination)) and (float(bin_info[-3]) > float(params.contamination)):
                    filter_h.write(info)
                    bin_path = os.path.join(params.fa_dir, bin_info[0] + ".fa")
                    cmd = "ln -s %s %s" % (bin_path, params.link_dir)
                    subprocess.run(cmd, shell=true)
'''
