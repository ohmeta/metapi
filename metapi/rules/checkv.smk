if config["params"]["checkv"]["do"]:
    rule checkv_download_db:
        output:
            expand([
                os.path.join(config["params"]["checkv"]["db"], "genome_db/checkv_{gfile}"),
                os.path.join(config["params"]["checkv"]["db"], "hmm_db/{hfile}")],
                gfile=["error.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.tsv", "source.tsv"],
                hfile=["checkv_hmms.tsv", "genome_lengths.tsv"])
        params:
            db = config["params"]["checkv"]["db"]
        benchmark:
            os.path.join(config["output"]["check"], "benchmark/checkv/checkv_download_db.benchmark.txt")
        log:
            os.path.join(config["output"]["check"], "logs/checkv/checkv_download_db.log")
        conda:
            config["envs"]["checkv"]
        shell:
            '''
            checkv download_database $(dirname {params.db}) \
            > {log} 2>&1
            '''
    

    rule checkv:
        input:
            db = expand([
                os.path.join(config["params"]["checkv"]["db"], "genome_db/checkv_{gfile}"),
                os.path.join(config["params"]["checkv"]["db"], "hmm_db/{hfile}")],
                gfile=["error.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.tsv", "source.tsv"],
                hfile=["checkv_hmms.tsv", "genome_lengths.tsv"]),
            viral = os.path.join(
                config["output"]["identify"],
                "vmags/{binning_group}.{assembly_group}.{assembler}/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.combined.fa")
        output:
            os.path.join(config["output"]["check"],
                         "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/checkv_done")
        benchmark:
            os.path.join(config["output"]["check"],
                         "benchmark/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.checkv.benchmark.txt")
        log:
            os.path.join(config["output"]["check"],
                         "logs/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.checkv.log")
        params:
            outdir = os.path.join(config["output"]["check"], "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}")
        conda:
            config["envs"]["checkv"]
        threads:
            config["params"]["checkv"]["threads"]
        shell:
            '''
            rm -rf {params.outdir}

            checkv end_to_end \
            {input.viral} \
            {params.outdir} \
            -t {threads} \
            -d {input.db} \
            >{log} 2>&1

            touch {output}
            '''

    #ASSEMBLY_GROUP = SAMPLES.reset_index().loc[:, ["assembly_group", "binning_group"]].drop_duplicates()

    #assembly_df_list = []
    #for assembler in ASSEMBLERS:
    #    assembly_df = ASSEMBLY_GROUP.copy()
    #    assembly_df["assembler"] = assembler
    #    assembly_df_list.append(assembly_df)
    #ASSEMBLY_GROUPS = pd.concat(assembly_df_list, axis=0)


    rule checkv_all:
        input:
            expand(
                os.path.join(config["output"]["check"],
                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/checkv_done"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"],
                identifier=config["params"]["checkv"]["checkv_identifier"])
 
else:
    rule checkv_all:
        input:


localrules:
    checkv_download_db