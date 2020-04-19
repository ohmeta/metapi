rule checkm_lineage_wf:
    input:
        bins_dir = os.path.join(
            config["output"]["binning"],
            "bins/{sample}.{assembler}.out/{binner}")
    output:
        table = os.path.join(
            config["output"]["checkm"],
            "table/{sample}/{sample}.{assembler}.{binner}.checkm.table.tsv"),
        data = os.path.join(
            config["output"]["checkm"],
            "data/{sample}/{sample}.{assembler}.{binner}.checkm.data.tar.gz")
    params:
        table_dir = os.path.join(config["output"]["checkm"], "table/{sample}"),
        data_dir = os.path.join(config["output"]["checkm"], "data/{sample}"),
        data_dir_ = os.path.join(config["output"]["checkm"],
                                 "data/{sample}/{sample}.{assembler}.{binner}")
    log:
        os.path.join(config["output"]["checkm"],
                     "logs/{sample}.{assembler}.{binner}.checkm.log")
    threads:
        config["params"]["checkm"]["threads"]
    shell:
        '''
        mkdir -p {params.table_dir}
        mkdir -p {params.data_dir}

        num=$(find {input.bins_dir} -type f -name "*.fa" | wc -l)

        if [[ $num > 0 ]]; then
            checkm lineage_wf \
            --tab_table \
            --file {output.table} \
            --threads {threads} \
            --extension fa \
            {input.bins_dir}/ \
            {params.data_dir_}/ > {log}
        else
            touch {output.table}
            mkdir -p {params.data_dir_}
        fi

        tar -czvf {output.data} {params.data_dir_}/
        rm -rf {params.data_dir_}
        '''


rule checkm_report:
    input:
        expand(
            os.path.join(
                config["output"]["checkm"],
                "table/{sample}/{sample}.{{assembler}}.{{binner}}.checkm.table.tsv"),
            sample=SAMPLES.index.unique())
    output:
        table = os.path.join(config["output"]["checkm"],
                             "report/{assembler}.{binner}.checkm.table.tsv")
    threads:
        config["params"]["checkm"]["threads"]
    run:
        metapi.checkm_report(input, output.table, threads)
       

rule checkm_link_bins:
    input:
        table = os.path.join(config["output"]["checkm"],
                             "report/{assembler}.{binner}.checkm.table.tsv")
    output:
        bins_hq_dir = directory(
            os.path.join(config["output"]["checkm"],
                         "bins_quality_high/{assembler}.{binner}.links")),
        bins_mq_dir = directory(
            os.path.join(config["output"]["checkm"],
                         "bins_quality_medium/{assembler}.{binner}.links")),
        bins_lq_dir = directory(
            os.path.join(config["output"]["checkm"],
                         "bins_quality_low/{assembler}.{binner}.links")),
        bins_hmq_dir = directory(
            os.path.join(config["output"]["checkm"],
                         "bins_quality_high_medium/{assembler}.{binner}.links"))
    params:
        bins_dir = os.path.join(config["output"]["binning"], "bins"),
        bin_suffix = ".fa",
        standard = config["params"]["checkm"]["standard"] + "_quality_level",
        assembler = "{assembler}",
        binner = "{binner}"
    run:
        if os.path.exists(output.bins_hq_dir):
            os.rmdir(output.bins_hq_dir)
        if os.path.exists(output.bins_mq_dir):
            os.rmdir(output.bins_mq_dir)
        if os.path.exists(output.bins_lq_dir):
            os.rmdir(output.bins_lq_dir)
        if os.path.exists(output.bins_hmq_dir):
            os.rmdir(output.bins_hmq_dir)

        os.mkdir(output.bins_hq_dir)
        os.mkdir(output.bins_mq_dir)
        os.mkdir(output.bins_lq_dir)
        os.mkdir(output.bins_hmq_dir)

        df = pd.read_csv(input.table, sep='\t').set_index("Bin Id")

        for bin_id in df.index:
            sample_id = bin_id.split(".")[0]
            bin_fa_path = os.path.realpath(
                os.path.join(
                    params.bins_dir,
                    sample_id + "." + params.assembler + ".out/" + \
                    params.binner + "/" + \
                    bin_id + params.bin_suffix))

            if df.loc[bin_id, params.standard] == "high_quality":
                os.symlink(bin_fa_path,
                           os.path.join(output.bins_hq_dir,
                                        bin_id + params.bin_suffix))
                os.symlink(bin_fa_path,
                           os.path.join(output.bins_hmq_dir,
                                        bin_id + params.bin_suffix))

            if df.loc[bin_id, params.standard] == "medium_quality":
                os.symlink(bin_fa_path,
                           os.path.join(output.bins_mq_dir,
                                        bin_id + params.bin_suffix))
                os.symlink(bin_fa_path,
                           os.path.join(output.bins_hmq_dir,
                                        bin_id + params.bin_suffix))

            if df.loc[bin_id, params.standard] == "low_quality":
                os.symlink(bin_fa_path,
                           os.path.join(output.bins_lq_dir,
                                        bin_id + params.bin_suffix))


rule checkm_all:
    input:
        expand([
            os.path.join(
                config["output"]["checkm"],
                "table/{sample}/{sample}.{assembler}.{binner}.checkm.table.tsv"),
            os.path.join(
                config["output"]["checkm"],
                "data/{sample}/{sample}.{assembler}.{binner}.checkm.data.tar.gz"),
            os.path.join(config["output"]["checkm"],
                         "report/{assembler}.{binner}.checkm.table.tsv"),
            os.path.join(config["output"]["checkm"],
                         "bins_quality_high/{assembler}.{binner}.links"),
            os.path.join(config["output"]["checkm"],
                         "bins_quality_medium/{assembler}.{binner}.links"),
            os.path.join(config["output"]["checkm"],
                         "bins_quality_low/{assembler}.{binner}.links"),
            os.path.join(config["output"]["checkm"],
                         "bins_quality_high_medium/{assembler}.{binner}.links")],
               assembler=ASSEMBLERS,
               binner=BINNERS,
               sample=SAMPLES.index.unique())
