checkpoint checkm_prepare:
    input:
        gene_table = os.path.join(
            config["output"]["predict"],
            "report/mags_gene_stats_{assembler}_{binner_checkm}.tsv.gz")
    output:
        mags_dir = directory(os.path.join(
            config["output"]["check"],
            "mags_input/{assembler}.{binner_checkm}"))
    params:
        batch_num = config["params"]["checkm"]["batch_num"]
    run:
        shell(f'''rm -rf {output.mags_dir}''')
        metapi.checkm_prepare(input.gene_table, params.batch_num, output.mags_dir)
        shell(f'''mkdir -p {output.mags_dir}/../../table/checkm''')
        shell(f'''mkdir -p {output.mags_dir}/../../data/checkm''')


rule checkm_lineage_wf:
    input:
        os.path.join(
            config["output"]["check"],
            "mags_input/{assembler}.{binner_checkm}/mags_input.{batchid}.tsv")
    output:
        table = os.path.join(
            config["output"]["check"],
            "table/checkm/checkm.table.{assembler}.{binner_checkm}.{batchid}.tsv.gz"),
        data = os.path.join(
            config["output"]["check"],
            "data/checkm/checkm.data.{assembler}.{binner_checkm}.{batchid}.tar.gz")
    log:
        os.path.join(
            config["output"]["check"],
            "logs/checkm_lineage_wf/{assembler}.{binner_checkm}.{batchid}.log")
    benchmark:
        os.path.join(
            config["output"]["check"],
            "benchmark/checkm_lineage_wf/{assembler}.{binner_checkm}.{batchid}.txt")
    wildcard_constraints:
        batchid="\d+"
    params:
        pplacer_threads = config["params"]["checkm"]["pplacer_threads"],
        reduced_tree = "--reduced_tree" if config["params"]["checkm"]["reduced_tree"] else "",
    threads:
        config["params"]["checkm"]["threads"]
    conda:
        config["envs"]["checkm"]
    shell:
        '''
        DATA={output.data}
        DATADIR=${{DATA%.tar.gz}}
        TABLEGZ={output.table}
        TABLE=${{TABLEGZ%.gz}}

        rm -rf $DATA $DATADIR $TABLEGZ $TABLE

        if [[ `wc -l {input} | awk '{{print $1}}'` -eq 0 ]];
        then
            echo "No genome found, please check the input again" > {log} 2>&1
            echo "Touch empty file and directory" >> {log} 2>&1

            touch $TABLE
            pigz -f $TABLE
            mkdir -p $DATADIR
            tar -czvf $DATA $DATADIR
            rm -rf $DATADIR
        else
            checkm lineage_wf \
            --tab_table \
            --file $TABLE \
            --threads {threads} \
            --pplacer_threads {params.pplacer_threads} \
            {params.reduced_tree} \
            --extension faa \
            --genes \
            {input} \
            $DATADIR \
            >{log} 2>&1

            pigz -f $TABLE
            tar -czvf $DATA $DATADIR
            rm -rf $DATADIR
        fi
        '''


def aggregate_checkm_output(wildcards):
    checkpoint_output = checkpoints.checkm_prepare.get(**wildcards).output.mags_dir

    return expand(os.path.join(
        config["output"]["check"],
        "table/checkm/checkm.table.{assembler}.{binner_checkm}.{batchid}.tsv.gz"),
        assembler=wildcards.assembler,
        binner_checkm=wildcards.binner_checkm,
        batchid=list(set([i.split("/")[0] \
        for i in glob_wildcards(os.path.join(
            checkpoint_output,
            "mags_input.{batchid}.tsv")).batchid])))


rule checkm_report:
    input:
        checkm_table = aggregate_checkm_output,
        gene_table = os.path.join(
            config["output"]["predict"],
            "report/mags_gene_stats_{assembler}_{binner_checkm}.tsv.gz"),
        mags_report = os.path.join(
            config["output"]["binning"],
            "report/assembly_stats_{assembler}_{binner_checkm}.tsv.gz")
    output:
        genomes_info = os.path.join(
            config["output"]["check"],
            "report/checkm/checkm_table_{assembler}_{binner_checkm}.tsv.gz"),
        mags_hq = os.path.join(
            config["output"]["check"],
            "report/checkm/MAGs_hq_{assembler}_{binner_checkm}.tsv.gz"),
        mags_mq = os.path.join(
            config["output"]["check"],
            "report/checkm/MAGs_mq_{assembler}_{binner_checkm}.tsv.gz"),
        mags_lq = os.path.join(
            config["output"]["check"],
            "report/checkm/MAGs_lq_{assembler}_{binner_checkm}.tsv.gz"),
        mags_hmq = os.path.join(
            config["output"]["check"],
            "report/checkm/MAGs_hmq_{assembler}_{binner_checkm}.tsv.gz")
    params:
        standard = config["params"]["checkm"]["standard"] + "_quality_level"
    threads:
        config["params"]["checkm"]["threads"]
    run:
        import pandas as pd

        checkm_table = metapi.checkm_reporter(input.checkm_table, None, threads)
        gene_table = pd.read_csv(input.gene_table, sep="\t")
        mags_report = metapi.extract_mags_report(input.mags_report)

        genomes_info = pd.merge(mags_report, gene_table, how="inner", on=["bin_id", "bin_file"])\
                            .merge(checkm_table, how="inner", on="bin_id")

        genomes_info["genome"] = genomes_info["bin_id"] + ".fa.gz"
        genomes_info.to_csv(output.genomes_info, sep='\t', index=False)

        genomes_info.query('%s=="high_quality"' % params.standard)\
            .loc[:, "bin_file"].to_csv(output.mags_hq, sep='\t', index=False, header=False)

        genomes_info.query('%s=="high_quality" or %s=="medium_quality"' % (params.standard, params.standard))\
            .loc[:, "bin_file"].to_csv(output.mags_hmq, sep='\t', index=False, header=False)

        genomes_info.query('%s=="medium_quality"' % params.standard)\
            .loc[:, "bin_file"].to_csv(output.mags_mq, sep='\t', index=False, header=False)

        genomes_info.query('%s=="low_quality"' % params.standard)\
            .loc[:, "bin_file"].to_csv(output.mags_lq, sep='\t', index=False, header=False)


rule checkm_report_merge:
    input:
        genomes_info = expand(os.path.join(
            config["output"]["check"],
            "report/checkm/checkm_table_{{assembler}}_{binner_checkm}.tsv.gz"),
            binner_checkm=BINNERS_CHECKM),
        mags_hmq = expand(os.path.join(
            config["output"]["check"],
            "report/checkm/MAGs_hmq_{{assembler}}_{binner_checkm}.tsv.gz"),
            binner_checkm=BINNERS_CHECKM)
    output:
        genomes_info = os.path.join(
            config["output"]["check"],
            "report/checkm/checkm_table_genomes_info.{assembler}.all.tsv"),
        genomes_info_simple = os.path.join(
            config["output"]["check"],
            "report/checkm/checkm_table_genomes_info.{assembler}.all.simple.csv"),
        mags_hmq = os.path.join(
            config["output"]["check"],
            "report/checkm/MAGs_hmq.{assembler}.all.tsv")
    run:
        import os
        import pandas as pd

        mags_hmq_list = [pd.read_csv(i, sep="\t", header=None) for i in input.mags_hmq]
        pd.concat(mags_hmq_list, axis=0).to_csv(output.mags_hmq, header=False, sep="\t", index=False)

        genomes_info_list = [pd.read_csv(i, sep="\t") for i in input.genomes_info]
        genomes_info_df = pd.concat(genomes_info_list)
        genomes_info_df.to_csv(output.genomes_info, sep="\t", index=False)

        genomes_info_df_simple = genomes_info_df.loc[:, ["bin_file", "completeness", "contamination"]]
        genomes_info_df_simple["genome"] = genomes_info_df_simple.apply(
            lambda x: os.path.splitext(x["bin_file"])[0], axis=1
        )
        genomes_info_df_simple\
            .loc[:, ["genome", "completeness", "contamination"]]\
            .to_csv(output.genomes_info_simple, index=False)



if config["params"]["checkm"]["do"]:
    rule checkm_all:
        input:
            expand([
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/checkm_table_{assembler}_{binner_checkm}.tsv.gz"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/MAGs_hq_{assembler}_{binner_checkm}.tsv.gz"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/MAGs_mq_{assembler}_{binner_checkm}.tsv.gz"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/MAGs_lq_{assembler}_{binner_checkm}.tsv.gz"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/MAGs_hmq_{assembler}_{binner_checkm}.tsv.gz"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/checkm_table_genomes_info.{assembler}.all.tsv"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/checkm_table_genomes_info.{assembler}.all.simple.csv"),
                os.path.join(
                    config["output"]["check"],
                    "report/checkm/MAGs_hmq.{assembler}.all.tsv")],
                    assembler=ASSEMBLERS,
                    binner_checkm=BINNERS_CHECKM)

else:
    rule checkm_all:
        input:


localrules:
    checkm_report_merge,
    checkm_prepare,
    checkm_report,
    checkm_all