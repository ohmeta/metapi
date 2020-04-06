rule checkm_lineage_wf:
    input:
        bins_dir = os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.{binner}_out")
    output:
        checkm_txt = os.path.join(config["results"]["checkm"]["out"], "{sample}.{assembler}.{binner}.checkm.txt"),
        checkm_targz = os.path.join(config["results"]["checkm"]["data"], "{sample}.{assembler}.{binner}.checkm.data.tar.gz")
    params:
        txt_dir = directory(config["results"]["checkm"]["out"]),
        data_dir = directory(config["results"]["checkm"]["data"]),
        checkm_data_dir = directory(os.path.join(config["results"]["checkm"]["data"], "{sample}.{assembler}.{binner}"))
    log:
        os.path.join(config["logs"]["checkm"], "{sample}.{assembler}.{binner}.checkm.log")
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
        expand("{checkmout}/{sample}.{assembler}.{binner}.checkm.txt",
               checkmout=config["results"]["checkm"]["out"],
               sample=_samples.index.unique(),
               assembler=config["params"]["assembler"],
               binner=config["params"]["binning"]["binner"])
    output:
        os.path.join(config["results"]["checkm"]["base_dir"],
                     "{assembler}.{binner}.checkm.out.tsv")
    threads:
        config["params"]["checkm"]["threads"]
    run:
        from metapi import checkmer
        checkmer.report(input, output[0], threads)
       

rule checkm_link_bins:
    input:
        os.path.join(config["results"]["checkm"]["base_dir"],
                     "{assembler}.{binner}.checkm.out.tsv")
    output:
        bins_hq = directory(os.path.join(config["results"]["checkm"]["base_dir"],
                                         "bins.{assembler}.{binner}_out.hq")),
        bins_mq = directory(os.path.join(config["results"]["checkm"]["base_dir"],
                                         "bins.{assembler}.{binner}_out.mq")),
        bins_lq = directory(os.path.join(config["results"]["checkm"]["base_dir"],
                                         "bins.{assembler}.{binner}_out.lq"))
    params:
        standard = config["params"]["checkm"]["standard"] + "_quality_level",
        bin_suffix = ".fa",
        bin_dir = config["results"]["binning"]["bins"],
        assembler = "{assembler}",
        binner = "{binner}"
    run:
        import pandas as pd
        import os

        os.rmdir(output.bins_hq)
        os.rmdir(output.bins_mq)
        os.rmdir(output.bins_lq)

        os.mkdir(output.bins_hq)
        os.mkdir(output.bins_mq)
        os.mkdir(output.bins_lq)

        df = pd.read_csv(input[0], sep='\t').set_index("bin_id")

        for bin_id in df.index:
            sample_id = bin_id.split(".")[0]
            bin_fa = os.path.join(os.path.join(params.bin_dir,
                                               sample_id + "." + params.assembler + "." + params.binner + "_out"),
                                  bin_id + params.bin_suffix)

            if df.loc[bin_id, params.standard] == "high_quality":
                os.symlink("../../../" + bin_fa, output.bins_hq + "/")

            if df.loc[bin_id, params.standard] == "medium_quality":
                os.symlink("../../../" + bin_fa, output.bins_mq + "/")

            if df.loc[bin_id, params.standard] == "low_quality":
                os.symlink("../../../" + bin_fa, output.bins_lq + "/")


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
