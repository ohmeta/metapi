if config["params"]["classify"]["gtdbtk"]["do"]:
    checkpoint coclassify_hmq_bins_gtdbtk_prepare:
        input:
            bins_hmq = os.path.join(
                config["output"]["cocheckm"],
                "report/{assembler_co}_{binner_checkm}_bins_hmq.tsv")
        output:
            out_dir = directory(
                os.path.join(config["output"]["coclassify"],
                             "bins_hmq_input/{assembler_co}_{binner_checkm}"))
        params:
            batch_num = config["params"]["classify"]["gtdbtk"]["batch_num"]
        run:
            import pandas as pd
            import os

            df = pd.read_csv(input.bins_hmq, names=["path"])
            df["id"] = df.apply(lambda x: os.path.basename(x["path"]), axis=1)

            os.makedirs(output.out_dir, exist_ok=True)

            if len(df) > 0:
                for batch_id in range(0, len(df), params.batch_num):
                    df_ = df.iloc[batch_id:batch_id + params.batch_num][["path", "id"]]
                    df_.to_csv(os.path.join(output.out_dir, "bins_hmq_%d.tsv" % batch_id),
                               sep='\t', index=False, header=None)
            else:
                shell('''touch {output.out_dir}/bins_hmq_0.tsv''')


    rule coclassify_hmq_bins_gtdbtk:
        input:
            bins_hmq = os.path.join(
                config["output"]["coclassify"],
                "bins_hmq_input/{assembler_co}_{binner_checkm}/bins_hmq_{batchid}.tsv")
        output:
            done = os.path.join(
                config["output"]["coclassify"],
                "table/{assembler_co}.{binner_checkm}.gtdbtk.out.{batchid}/done")
        wildcard_constraints:
            batchid="\d+"
        log:
            os.path.join(
                config["output"]["coclassify"],
                "logs/bins_hmq_{batchid}.{assembler_co}.{binner_checkm}.gtdbtk.log")
        params:
            bin_suffix = config["params"]["binning"]["bin_suffix"],
            out_dir = os.path.join(
                config["output"]["coclassify"],
                "table/{assembler}.{binner_checkm}.gtdbtk.out.{batchid}"),
            pplacer_threads = config["params"]["classify"]["gtdbtk"]["pplacer_threads"]
        threads:
            config["params"]["classify"]["threads"]
        shell:
            '''
            gtdbtk classify_wf \
            --batchfile {input.bins_hmq} \
            --out_dir {params.out_dir} \
            --extension {params.bin_suffix} \
            --cpus {threads} \
            --pplacer_cpus {params.pplacer_threads} \
            > {log}

            touch {output.done}
            '''


    def aggregate_coclassify_gtdbtk_report_input(wildcards):
        checkpoint_output = checkpoints.coclassify_hmq_bins_gtdbtk_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["coclassify"],
            "table/{assembler_co}.{binner_checkm}.gtdbtk.out.{batchid}/done"),
                      assembler_co=wildcards.assembler_co,
                      binner_checkm=wildcards.binner_checkm,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(
                                                os.path.join(checkpoint_output,
                                                             "bins_hmq_{batchid}.tsv")).batchid])))


    rule coclassify_hmq_bins_gtdbtk_report:
        input:
            aggregate_coclassify_gtdbtk_report_input
        output:
            table_gtdb_a = os.path.join(
                config["output"]["coclassify"],
                "report/bins_hmq.{assembler_co}.{binner_checkm}.gtdbtk.archaea.gtdb.tsv"),
            table_gtdb_b = os.path.join(
                config["output"]["coclassify"],
                "report/bins_hmq.{assembler_co}.{binner_checkm}.gtdbtk.bacteria.gtdb.tsv")
        threads:
            8
        run:
            import os

            ar122_list = []
            bac120_list = []

            for i in input:
                ar122_tsv = os.path.join(os.path.dirname(i), "gtdbtk.ar122.summary.tsv")
                bac120_tsv = os.path.join(os.path.dirname(i), "gtdbtk.bac120.summary.tsv")
                if os.path.exists(ar122_tsv):
                    ar122_list.append(ar122_tsv)
                if os.path.exists(bac120_tsv):
                    bac120_list.append(bac120_tsv)

            metapi.merge(ar122_list, metapi.parse, threads, output=output.table_gtdb_a)
            metapi.merge(bac120_list, metapi.parse, threads, output=output.table_gtdb_b)


    rule coclassify_hmq_bins_gtdbtk_all:
        input:
            expand(
                os.path.join(
                    config["output"]["coclassify"],
                    "report/bins_hmq.{assembler_co}.{binner_checkm}.gtdbtk.{taxonomy}.{system}.tsv"),
                taxonomy=["archaea", "bacteria"],
                system=["gtdb"],
                assembler_co=ASSEMBLERS_CO,
                binner_checkm=BINNERS_CHECKM),

            rules.cocheckm_all.input,

else:
    rule coclassify_hmq_bins_gtdbtk_all:
        input:


rule coclassify_all:
    input:
        rules.coclassify_hmq_bins_gtdbtk_all.input


rule classify_hmq_bins_gtdbtk_all:
    input:
        rules.single_classify_hmq_bins_gtdbtk_all.input,
        rules.coclassify_hmq_bins_gtdbtk_all.input


rule classify_all:
    input:
        rules.single_classify_all.input,
        rules.coclassify_all.input
