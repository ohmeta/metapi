if config["params"]["checkm"]["do"]:
    checkpoint checkm_prepare:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembly_group}.{{assembler}}.prodigal.out/{{binner_checkm}}/predict_done"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
        output:
            bins_dir = directory(os.path.join(config["output"]["checkm"],
                                               "bins_input/{assembler}.{binner_checkm}.bins_input"))
        params:
            batch_num = config["params"]["checkm"]["batch_num"]
        run:
            metapi.checkm_prepare(input, params.batch_num, output.bins_dir)


    rule checkm_lineage_wf:
        input:
            os.path.join(config["output"]["checkm"],
                         "bins_input/{assembler}.{binner_checkm}.bins_input/bins_input.{batchid}.tsv")
        output:
            table = os.path.join(config["output"]["checkm"],
                                 "table/checkm.table.{batchid}.{assembler}.{binner_checkm}.tsv"),
            data = directory(os.path.join(config["output"]["checkm"],
                                          "data/checkm.data.{batchid}.{assembler}.{binner_checkm}"))
        wildcard_constraints:
            batchid="\d+"
        params:
            pplacer_threads = config["params"]["checkm"]["pplacer_threads"],
            reduced_tree = "--reduced_tree" if config["params"]["checkm"]["reduced_tree"] else "",
        conda:
            config["envs"]["checkm"]
        log:
            os.path.join(config["output"]["checkm"],
                "logs/{assembler}.{binner_checkm}.checkm.{batchid}.log")
        benchmark:
            os.path.join(config["output"]["checkm"],
                         "benchmark/{assembler}.{binner_checkm}.checkm.{batchid}.benchmark.txt")
        threads:
            config["params"]["checkm"]["threads"]
        shell:
            '''
            checkm lineage_wf \
            --tab_table \
            --file {output.table} \
            --threads {threads} \
            --pplacer_threads {params.pplacer_threads} \
            {params.reduced_tree} \
            --extension faa \
            --genes \
            {input} \
            {output.data} \
            > {log}
            '''


    def aggregate_checkm_output(wildcards):
        checkpoint_output = checkpoints.checkm_prepare.get(**wildcards).output.bins_dir

        return expand(os.path.join(
            config["output"]["checkm"],
            "table/checkm.table.{batchid}.{assembler}.{binner_checkm}.tsv"),
                      assembler=wildcards.assembler,
                      binner_checkm=wildcards.binner_checkm,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(os.path.join(checkpoint_output,
                                                                             "bins_input.{batchid}.tsv")).batchid])))

   
    rule checkm_report:
        input:
            aggregate_checkm_output
        output:
            table = os.path.join(config["output"]["checkm"],
                                 "report/checkm_table_{assembler}_{binner_checkm}.tsv"),
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
            standard = config["params"]["checkm"]["standard"] + "_quality_level",
            assembler = "{assembler}",
            binner = "{binner_checkm}"
        run:
            import pandas as pd

            df = metapi.checkm_reporter(input, output.table, threads)

            def get_bin_path(row):
                bin_fa_path = os.path.realpath(
                    os.path.join(
                        params.bins_dir,
                        "%s.%s.out/%s/%s.fa" % (
                            row["bin_id"].split(".")[0],
                            params.assembler,
                            params.binner,
                            row["bin_id"])))
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


    rule checkm_all:
        input:
            expand([
                os.path.join(config["output"]["checkm"],
                             "report/checkm_table_{assembler}_{binner_checkm}.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_hq.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_mq.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_lq.tsv"),
                os.path.join(config["output"]["checkm"],
                             "report/{assembler}_{binner_checkm}_bins_hmq.tsv")],
                assembler=ASSEMBLERS,
                binner_checkm=BINNERS_CHECKM)

            #rules.predict_bins_gene_prodigal_all.input,
            #rules.binning_all.input,

else:
    rule checkm_all:
        input: