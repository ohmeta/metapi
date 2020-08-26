if config["params"]["checkm"]["do"]:
    checkpoint checkm_prepare:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "bins_gene/{{assembler}}.{{binner_checkm}}.prodigal.out/{sample}/done"),
                   sample=SAMPLES.index.unique())
        output:
            directory(os.path.join(
                config["output"]["checkm"],
                "bins_input/{assembler}.{binner_checkm}.links"))
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


    rule checkm_lineage_wf:
        input:
            os.path.join(config["output"]["checkm"],
                         "bins_input/{assembler}.{binner_checkm}.links/bins_{batchid}")
        output:
            table = os.path.join(
                config["output"]["checkm"],
                "table/bins_{batchid}/bins_{batchid}.{assembler}.{binner_checkm}.checkm.table.tsv"),
            data = os.path.join(
                config["output"]["checkm"],
                "data/bins_{batchid}/bins_{batchid}.{assembler}.{binner_checkm}.checkm.data.tar.gz")
        wildcard_constraints:
            batchid="\d+"
        params:
            suffix = "faa",
            pplacer_threads = config["params"]["checkm"]["pplacer_threads"],
            reduced_tree = "--reduced_tree" if config["params"]["checkm"]["reduced_tree"] else "",
            table_dir = os.path.join(config["output"]["checkm"], "table/bins_{batchid}"),
            data_dir = os.path.join(config["output"]["checkm"], "data/bins_{batchid}"),
            data_dir_temp = os.path.join(
                config["output"]["checkm"],
                "data/bins_{batchid}/bins_{batchid}.{assembler}.{binner_checkm}")
        log:
            os.path.join(
                config["output"]["checkm"],
                "logs/bins_{batchid}.{assembler}.{binner_checkm}.checkm.log")
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
                    --pplacer_threads {params.pplacer_threads} \
                    {params.reduced_tree} \
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


    def aggregate_checkm_report_input(wildcards):
        checkpoint_output = checkpoints.checkm_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["checkm"],
            "table/bins_{batchid}/bins_{batchid}.{assembler}.{binner_checkm}.checkm.table.tsv"),
                      assembler=wildcards.assembler,
                      binner_checkm=wildcards.binner_checkm,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(os.path.join(checkpoint_output,
                                                                             "bins_{batchid}")).batchid])))

   
    rule checkm_report:
        input:
            aggregate_checkm_report_input
        output:
            table = os.path.join(config["output"]["checkm"],
                                 "report/{assembler}_{binner_checkm}_checkm_table.tsv"),
            bins_hq = os.path.join(config["output"]["checkm"],
                                   "report/{assembler}_{binner_checkm}_bins_hq.tsv"),
            bins_mq = os.path.join(config["output"]["checkm"],
                                   "report/{assembler}_{binner_checkm}_bins_mq.tsv"),
            bins_lq = os.path.join(config["output"]["checkm"],
                                   "report/{assembler}_{binner_checkm}_bins_lq.tsv"),
            bins_hmq = os.path.join(config["output"]["checkm"],
                                    "report/{assembler}_{binner_checkm}_bins_hmq.tsv")
        threads:
            config["params"]["checkm"]["threads"]
        params:
            bins_dir = os.path.join(config["output"]["binning"], "bins"),
            bin_suffix = config["params"]["binning"]["bin_suffix"],
            standard = config["params"]["checkm"]["standard"] + "_quality_level",
            assembler = "{assembler}",
            binner = "{binner_checkm}"
        run:
            import pandas as pd

            df = metapi.checkm_report(input, output.table, threads)

            def get_bin_path(row):
                bin_fa_path = os.path.realpath(
                    os.path.join(
                        params.bins_dir,
                        "%s.%s.out/%s/%s.%s" % (
                            row["bin_id"].split(".")[0],
                            params.assembler,
                            params.binner,
                            row["bin_id"],
                            params.bin_suffix)))
                return bin_fa_path

            df["bin_fa_path"] = df.apply(lambda x: get_bin_path(x), axis=1)

            df.query('%s=="high_quality"' % params.standard)\
              .loc[:, "bin_fa_path"]\
              .to_csv(output.bins_hq, sep='\t', index=False, header=False)

            df.query('%s=="high_quality" or %s=="medium_quality"' % (params.standard, params.standard))\
              .loc[:, "bin_fa_path"]\
              .to_csv(output.bins_hmq, sep='\t', index=False, header=False)

            df.query('%s=="medium_quality"' % params.standard)\
              .loc[:, "bin_fa_path"]\
              .to_csv(output.bins_mq, sep='\t', index=False, header=False)

            df.query('%s=="low_quality"' % params.standard)\
              .loc[:, "bin_fa_path"]\
              .to_csv(output.bins_lq, sep='\t', index=False, header=False)


    rule single_checkm_all:
        input:
            expand([
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_checkm_table.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_hq.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_mq.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_lq.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_hmq.tsv")],
                assembler=ASSEMBLERS,
                binner_checkm=BINNERS_CHECKM),

            rules.predict_bins_gene_prodigal_all.input,
            rules.binning_all.input,

else:
    rule single_checkm_all:
        input:
