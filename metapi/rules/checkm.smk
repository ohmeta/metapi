if config["params"]["checkm"]["do"]:
    checkpoint checkm_prepare:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "bins_gene/{{assembler}}.{{binner}}.prodigal.out/{sample}/done"),
                   sample=SAMPLES.index.unique()) \
                if config["params"]["predict"]["bins_to_gene"]["prodigal"]["do"] else \
                   expand(os.path.join(
                       config["output"]["binning"],
                       "bins/{sample}.{{assembler}}.out/{{binner}}"),
                          sample=SAMPLES.index.unique())
        output:
            directory(os.path.join(
                config["output"]["checkm"],
                "bins_input/{assembler}.{binner}.links"))
        params:
            suffix = "faa" \
                if config["params"]["predict"]["bins_to_gene"]["prodigal"]["do"] \
                   else "fa",
             batch_num = config["params"]["checkm"]["batch_num"]
        run:
            import os
            import glob
            import pprint

            if os.path.exists(output[0]):
                os.rmdir(output[0])

            bin_list = []
            if params.suffix == "faa":
                for i in input:
                   bin_list += [os.path.realpath(j) \
                                for j in glob.glob(os.path.join(os.path.dirname(i), "*.faa"))]
            if params.suffix == "fa":
                for i in input:
                    bin_list += [os.path.realpath(j) \
                                 for j in glob.glob(os.path.join(i, "*.fa"))]

            if len(bin_list) > 0:
                for batch_id in range(0, len(bin_list), params.batch_num):
                    batch_dir = os.path.join(output[0], "bins_%d" % batch_id)
                    os.makedirs(batch_dir, exist_ok=True)

                    for bin_file in bin_list[batch_id:batch_id + params.batch_num]:
                        os.symlink(bin_file,
                                   os.path.join(batch_dir,
                                                os.path.basename(bin_file)))
            else:
                os.makedirs(os.path.join(output[0], "bins_0"), exist_ok=True)


    rule checkm_lineage_wf:
        input:
            os.path.join(config["output"]["checkm"],
                         "bins_input/{assembler}.{binner}.links/bins_{batchid}")
        output:
            table = os.path.join(
                config["output"]["checkm"],
                "table/bins_{batchid}/bins_{batchid}.{assembler}.{binner}.checkm.table.tsv"),
            data = os.path.join(
                config["output"]["checkm"],
                "data/bins_{batchid}/bins_{batchid}.{assembler}.{binner}.checkm.data.tar.gz")
        wildcard_constraints:
            batchid="\d+"
        params:
            suffix = "faa" \
                if config["params"]["predict"]["bins_to_gene"]["prodigal"]["do"] \
                   else "fa",
            genes = "--genes" \
                if config["params"]["predict"]["bins_to_gene"]["prodigal"]["do"] \
                   else "",
            table_dir = os.path.join(config["output"]["checkm"], "table/bins_{batchid}"),
            data_dir = os.path.join(config["output"]["checkm"], "data/bins_{batchid}"),
            data_dir_temp = os.path.join(config["output"]["checkm"],
                                     "data/bins_{batchid}/bins_{batchid}.{assembler}.{binner}")
        log:
            os.path.join(
                config["output"]["checkm"],
                "logs/bins_{batchid}.{assembler}.{binner}.checkm.log")
        threads:
            config["params"]["checkm"]["threads"]
        run:
            import os
            import glob

            os.makedirs(params.table_dir, exist_ok=True)
            os.makedirs(params.data_dir, exist_ok=True)

            count = len(glob.glob(os.path.join(input[0], "*.%s" % params.suffix)))

            if count > 0:
                shell(
                    '''
                    checkm lineage_wf \
                    --tab_table \
                    --file {output.table} \
                    --threads {threads} \
                    --extension {params.suffix} \
                    {params.genes} \
                    {input}/ \
                    {params.data_dir_temp}/ > {log}
                    ''')
            else:
                shell('''touch {output.table}''')
                shell('''mkdir -p {params.data_dir_temp}''')

            shell('''tar -czvf {output.data} {params.data_dir_temp}/''')
            shell('''rm -rf {params.data_dir_temp}''')


    def aggregate_checkm_report_input(wildcards):
        checkpoint_output = checkpoints.checkm_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["checkm"],
            "table/bins_{batchid}/bins_{batchid}.{assembler}.{binner}.checkm.table.tsv"),
                      assembler=wildcards.assembler,
                      binner=wildcards.binner,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(os.path.join(checkpoint_output,
                                                                             "bins_{batchid}")).batchid])))

   
    rule checkm_report:
        input:
            aggregate_checkm_report_input
        output:
            table = os.path.join(config["output"]["checkm"],
                                 "report/{assembler}_{binner}_checkm_table.tsv")
        threads:
            config["params"]["checkm"]["threads"]
        run:
            metapi.checkm_report(input, output.table, threads)
       

    rule checkm_link_bins:
        input:
            table = os.path.join(config["output"]["checkm"],
                                 "report/{assembler}_{binner}_checkm_table.tsv")
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
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner}_checkm_table.tsv"),
                os.path.join(config["output"]["checkm"],
                             "bins_hq/{assembler}.{binner}.links"),
                os.path.join(config["output"]["checkm"],
                             "bins_mq/{assembler}.{binner}.links"),
                os.path.join(config["output"]["checkm"],
                             "bins_lq/{assembler}.{binner}.links"),
                os.path.join(config["output"]["checkm"],
                             "bins_hmq/{assembler}.{binner}.links")],
                   assembler=ASSEMBLERS,
                   binner=BINNERS),

            rules.binning_all.input

else:
    rule checkm_all:
        input:
