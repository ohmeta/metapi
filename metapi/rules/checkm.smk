rule checkm_lineage_wf:
    input:
        bins_dir = os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.metabat2_out")
    output:
        checkm_txt = os.path.join(config["results"]["checkm"]["out"], "{sample}.{assembler}.checkm.txt"),
        checkm_targz = os.path.join(config["results"]["checkm"]["data"], "{sample}.{assembler}.checkm.data.tar.gz")
    params:
        txt_dir = directory(config["results"]["checkm"]["out"]),
        data_dir = directory(config["results"]["checkm"]["data"]),
        checkm_data_dir = directory(os.path.join(config["results"]["checkm"]["data"], "{sample}.{assembler}"))
    log:
        os.path.join(config["logs"]["checkm"], "{sample}.{assembler}.checkm.log")
    threads:
        config["params"]["checkm"]["threads"]
    shell:
        '''
        mkdir -p {params.txt_dir}
        mkdir -p {params.data_dir}
        mkdir -p {params.checkm_data_dir}
        num=$(find {input.bins_dir} -type f -name "*.fa" | wc -l)
        if [[ $num > 0 ]]; then
            rm -rf {params.checkm_data_dir}
            checkm lineage_wf -f {output.checkm_txt} -t {threads} -x fa {input.bins_dir}/ {params.checkm_data_dir}/ 2> {log}
        else
            echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" > {output.checkm_txt}
            echo "Bin Id                                              Marker lineage             # genomes   # markers   # marker sets    0     1     2    3    4    5+   Completeness   Contamination   Strain heterogeneity    " >> {output.checkm_txt} 
            echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> {output.checkm_txt}
        fi
        tar -czvf {output.checkm_targz} {params.checkm_data_dir}/
        rm -rf {params.checkm_data_dir}
        '''


rule checkm_report:
    input:
        expand("{checkmout}/{sample}.{assembler}.checkm.txt",
               checkmout=config["results"]["checkm"]["out"],
               sample=_samples.index.unique(),
               assembler=config["params"]["assembler"])
    output:
        os.path.join(config["results"]["checkm"]["basedir"], "checkm_out.tsv")
    params:
        completeness = config["params"]["checkm"]["completeness"],
        contamination = config["params"]["checkm"]["contamination"]
    threads:
        config["params"]["checkm"]["threads"]
    run:
        from metapi import checkmer
        checkmer.report(input, output[0], params.completeness, params.contamination)
       

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
