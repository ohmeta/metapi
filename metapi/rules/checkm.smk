rule checkm_lineage_wf:
    input:
        bins_dir = os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.metabat2_out")
    output:
        checkm_txt = os.path.join(config["results"]["checkm"]["out"], "{sample}.{assembler}.checkm.txt"),
        checkm_data_dir = directory(os.path.join(config["results"]["checkm"]["data"], "{sample}.{assembler}"))
    params:
        txt_dir = directory(config["results"]["checkm"]["out"]),
        data_dir = directory(config["results"]["checkm"]["data"]),
        checkm_env = config["params"]["checkm"]["env"]
    threads:
        config["params"]["checkm"]["threads"]
    shell:
        '''
        mkdir -p {params.txt_dir}
        mkdir -p {params.data_dir}
        num=$(find {input.bins_dir} -type f -name "*.fa" | wc -l)
        if [[ $num > 0 ]]; then
            set +u; source activate {params.checkm_env}; set -u;    
            checkm lineage_wf -f {output.checkm_txt} -t {threads} -x fa {input.bins_dir}/ {output.checkm_data_dir}/
        else
            echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" > {output.checkm_txt}
            echo "Bin Id                                              Marker lineage             # genomes   # markers   # marker sets    0     1     2    3    4    5+   Completeness   Contamination   Strain heterogeneity    " >> {output.checkm_txt} 
            echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> {output.checkm_txt}
        fi
        '''


# rule checkm_coverage:
#     input:
#         bins_dir = os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.metabat2_out"),
#         bam = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam"),
#         bai = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam.bai")
#     output:
#         os.path.join(config["results"]["checkm"]["coverage"], "{sample}.{assembler}.checkm_coverage.tsv")
#     log:
#         os.path.join(config["logs"]["checkm"], "{sample}.{assembler}.checkm_coverage.log")
#     params:
#         checkm_env = config["params"]["checkm"]["env"]
#     threads:
#         config["params"]["checkm"]["threads"]
#     shell:
#         '''
#         set +u; source activate {params.checkm_env}; set -u;
#         checkm coverage -x fa {input.bins_dir} {output} {input.bam} -t {threads} 2> {log}
#         '''

# rule checkm_profile:
#     input:
#         os.path.join(config["results"]["checkm"]["coverage"], "{sample}.{assembler}.checkm_coverage.tsv")
#     output:
#         os.path.join(config["results"]["checkm"]["profile"], "{sample}.{assembler}.checkm_profile.tsv")
#     log:
#         os.path.join(config["logs"]["checkm"], "{sample}.{assembler}.checkm_profile.log")
#     params:
#         checkm_env = config["params"]["checkm"]["env"]
#     shell:
#         '''
#         set +u; source activate {params.checkm_env}; set -u;
#         checkm profile -f {output} {input} 2> {log}
#         '''