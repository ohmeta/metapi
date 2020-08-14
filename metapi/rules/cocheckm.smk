if config["params"]["checkm"]["do"]:
    checkpoint cocheckm_prepare:
        input:
            os.path.join(
                config["output"]["copredict"],
                "bins_gene/{assembler_co}.{binner_checkm}.prodigal.out/all/done")
        output:
            directory(os.path.join(
                config["output"]["cocheckm"],
                "bins_input/{assembler_co}.{binner_checkm}.links"))
        params:
            suffix = "faa",
            batch_num = config["params"]["checkm"]["batch_num"]
        run:
            import os
            import glob
            import pprint

            if os.path.exists(output[0]):
                os.rmdir(output[0])

            bin_list = []
            for i in input:
                bin_list += [os.path.realpath(j) \
                             for j in glob.glob(os.path.join(os.path.dirname(i), "*.faa"))]

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


    rule cocheckm_lineage_wf:
        input:
            os.path.join(config["output"]["cocheckm"],
                         "bins_input/{assembler_co}.{binner_checkm}.links/bins_{batchid}")
        output:
            table = os.path.join(
                config["output"]["cocheckm"],
                "table/bins_{batchid}/bins_{batchid}.{assembler_co}.{binner_checkm}.checkm.table.tsv"),
            data = os.path.join(
                config["output"]["cocheckm"],
                "data/bins_{batchid}/bins_{batchid}.{assembler_co}.{binner_checkm}.checkm.data.tar.gz")
        wildcard_constraints:
            batchid="\d+"
        params:
            suffix = "faa",
            table_dir = os.path.join(config["output"]["cocheckm"], "table/bins_{batchid}"),
            data_dir = os.path.join(config["output"]["cocheckm"], "data/bins_{batchid}"),
            data_dir_temp = os.path.join(
                config["output"]["cocheckm"],
                "data/bins_{batchid}/bins_{batchid}.{assembler_co}.{binner_checkm}")
        log:
            os.path.join(
                config["output"]["cocheckm"],
                "logs/bins_{batchid}.{assembler_co}.{binner_checkm}.checkm.log")
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
                    --genes \
                    {input}/ \
                    {params.data_dir_temp}/ > {log}
                    ''')
            else:
                shell('''touch {output.table}''')
                shell('''mkdir -p {params.data_dir_temp}''')

            shell('''tar -czvf {output.data} {params.data_dir_temp}/''')
            shell('''rm -rf {params.data_dir_temp}''')


    def aggregate_cocheckm_report_input(wildcards):
        checkpoint_output = checkpoints.cocheckm_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["cocheckm"],
            "table/bins_{batchid}/bins_{batchid}.{assembler_co}.{binner_checkm}.checkm.table.tsv"),
                      assembler_co=wildcards.assembler_co,
                      binner_checkm=wildcards.binner_checkm,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(os.path.join(checkpoint_output,
                                                                             "bins_{batchid}")).batchid])))


    rule cocheckm_report:
        input:
            aggregate_cocheckm_report_input
        output:
            table = os.path.join(config["output"]["cocheckm"],
                                 "report/{assembler_co}_{binner_checkm}_checkm_table.tsv")
        threads:
            config["params"]["checkm"]["threads"]
        run:
            metapi.checkm_report(input, output.table, threads)


    rule cocheckm_link_bins:
        input:
            table = os.path.join(config["output"]["cocheckm"],
                                 "report/{assembler_co}_{binner_checkm}_checkm_table.tsv")
        output:
            bins_dir_hq = directory(
                os.path.join(config["output"]["cocheckm"],
                             "bins_hq/{assembler_co}.{binner_checkm}.links")),
            bins_dir_mq = directory(
                os.path.join(config["output"]["cocheckm"],
                             "bins_mq/{assembler_co}.{binner_checkm}.links")),
            bins_dir_lq = directory(
                os.path.join(config["output"]["cocheckm"],
                             "bins_lq/{assembler_co}.{binner_checkm}.links")),
            bins_dir_hmq = directory(
                os.path.join(config["output"]["cocheckm"],
                             "bins_hmq/{assembler_co}.{binner_checkm}.links"))
        params:
            bins_dir = os.path.join(config["output"]["cobinning"], "bins"),
            bin_suffix = ".fa",
            standard = config["params"]["checkm"]["standard"] + "_quality_level",
            assembler_co = "{assembler_co}",
            binner = "{binner_checkm}"
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

            df = pd.read_csv(input.table, sep='\t').set_index("bin_id")

            for bin_id in df.index:
                sample_id = bin_id.split(".")[0]
                bin_fa_path = os.path.realpath(
                    os.path.join(
                        params.bins_dir,
                        sample_id + "." + params.assembler_co + ".out/" + \
                        params.binner + "/" + \
                        bin_id + params.bin_suffix))

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


    rule cocheckm_all:
        input:
            expand([
                os.path.join(config["output"]["cocheckm"],
                             "report/{assembler_co}_{binner_checkm}_checkm_table.tsv"),
                os.path.join(config["output"]["cocheckm"],
                             "bins_hq/{assembler_co}.{binner_checkm}.links"),
                os.path.join(config["output"]["cocheckm"],
                             "bins_mq/{assembler_co}.{binner_checkm}.links"),
                os.path.join(config["output"]["cocheckm"],
                             "bins_lq/{assembler_co}.{binner_checkm}.links"),
                os.path.join(config["output"]["cocheckm"],
                             "bins_hmq/{assembler_co}.{binner_checkm}.links")],
                   assembler_co=ASSEMBLERS_CO,
                   binner_checkm=BINNERS_CHECKM),

            rules.copredict_bins_gene_prodigal_all.input,
            rules.cobinning_all.input,

else:
    rule cocheckm_all:
        input:


rule checkm_all:
    input:
        rules.single_checkm_all.input,
        rules.cocheckm_all.input
