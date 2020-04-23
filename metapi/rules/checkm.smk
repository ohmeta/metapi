if config["params"]["checkm"]["do"]:
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
            bin_suffix = "fa",
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
                --extension {params.bin_suffix} \
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
            bins_dir_hq = directory(
                os.path.join(config["output"]["checkm"],
                             "bins_hq/{assembler}.{binner}.links")),
            bins_dir_mq = directory(
                os.path.join(config["output"]["checkm"],
                             "bins_mq/{assembler}.{binner}.links")),
            bins_dir_lq = directory(
                os.path.join(config["output"]["checkm"],
                             "bins_lq/{assembler}.{binner}.links")),
            bins_dir_hmq = directory(
                os.path.join(config["output"]["checkm"],
                             "bins_hmq/{assembler}.{binner}.links"))
        params:
            bins_dir = os.path.join(config["output"]["binning"], "bins"),
            bin_suffix = "fa",
            standard = config["params"]["checkm"]["standard"] + "_quality_level",
            assembler = "{assembler}",
            binner = "{binner}"
        run:
            if os.path.exists(output.bins_dir_hq):
                os.rmdir(output.bins_dir_hq)
            if os.path.exists(output.bins_dir_mq):
                os.rmdir(output.bins_dir_mq)
            if os.path.exists(output.bins_dir_lq):
                os.rmdir(output.bins_dir_lq)
            if os.path.exists(output.bins_dir_hmq):
                os.rmdir(output.bins_dir_hmq)

            os.mkdir(output.bins_dir_hq)
            os.mkdir(output.bins_dir_mq)
            os.mkdir(output.bins_dir_lq)
            os.mkdir(output.bins_dir_hmq)

            df = pd.read_csv(input.table, sep='\t').set_index("Bin Id")

            for bin_id in df.index:
                sample_id = bin_id.split(".")[0]
                bin_fa_path = os.path.realpath(
                    os.path.join(
                        params.bins_dir,
                        sample_id + "." + params.assembler + ".out/" + \
                        params.binner + "/" + \
                        bin_id + "." + params.bin_suffix))

                if df.loc[bin_id, params.standard] == "high_quality":
                    os.symlink(bin_fa_path,
                               os.path.join(output.bins_dir_hq,
                                            bin_id + params.bin_suffix))
                    os.symlink(bin_fa_path,
                               os.path.join(output.bins_dir_hmq,
                                            bin_id + params.bin_suffix))

                if df.loc[bin_id, params.standard] == "medium_quality":
                    os.symlink(bin_fa_path,
                               os.path.join(output.bins_dir_mq,
                                            bin_id + params.bin_suffix))
                    os.symlink(bin_fa_path,
                               os.path.join(output.bins_dir_hmq,
                                            bin_id + params.bin_suffix))

                if df.loc[bin_id, params.standard] == "low_quality":
                    os.symlink(bin_fa_path,
                               os.path.join(output.bins_dir_lq,
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
                             "bins_hq/{assembler}.{binner}.links"),
                os.path.join(config["output"]["checkm"],
                             "bins_mq/{assembler}.{binner}.links"),
                os.path.join(config["output"]["checkm"],
                             "bins_lq/{assembler}.{binner}.links"),
                os.path.join(config["output"]["checkm"],
                             "bins_hmq/{assembler}.{binner}.links")],
                   assembler=ASSEMBLERS,
                   binner=BINNERS,
                   sample=SAMPLES.index.unique())

else:
    rule checkm_all:
        input:
