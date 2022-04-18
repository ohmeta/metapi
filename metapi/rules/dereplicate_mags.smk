rule dereplicate_mags_prepare:
    input:
        genomes_info = expand(os.path.join(config["output"]["checkm"],
                                           "report/checkm_table_{assembler}_{binner_checkm}.tsv"),
                              assembler=ASSEMBLERS,
                              binner_checkm=BINNERS_CHECKM),
        bins_hmq = expand(os.path.join(config["output"]["checkm"],
                                       "report/{assembler}_{binner_checkm}_bins_hmq.tsv"),
                          assembler=ASSEMBLERS,
                          binner_checkm=BINNERS_CHECKM)
    output:
        genomes_info = os.path.join(config["output"]["dereplicate"],
                                    "genomes_info/checkm_table_genomes_info.all.tsv"),
        bins_hmq = os.path.join(config["output"]["dereplicate"],
                                "genomes_info/bins_hmq.all.tsv")
    run:
        import pandas as pd

        genomes_info_list = [pd.read_csv(i, sep="\t") for i in input.genomes_info]
        pd.concat(genomes_info_list, axis=0).to_csv(output.genomes_info, sep="\t", index=False)

        bins_hmq_list = [pd.read_csv(i, sep="\t", header=None) for i in input.bins_hmq]
        pd.concat(bins_hmq_list, axis=0).to_csv(output.bins_hmq, header=False, sep="\t", index=False)


if config["params"]["dereplicate"]["drep"]["do"]:
    rule dereplicate_mags_drep:
        input:
            genomes_info = os.path.join(config["output"]["dereplicate"],
                                        "genomes_info/checkm_table_genomes_info.all.tsv"),
            bins_hmq = os.path.join(config["output"]["dereplicate"],
                                    "genomes_info/bins_hmq.all.tsv")
        output:
            os.path.join(config["output"]["dereplicate"], "genomes/hmq.bins.drep.out/drep_done")
        log:
            os.path.join(config["output"]["dereplicate"], "logs/hmq.bins.drep.log")
        benchmark:
            os.path.join(config["output"]["dereplicate"], "benchmark/drep.benchmark.txt")
        conda:
            config["envs"]["drep"]
        params:
            output_dir = os.path.join(config["output"]["dereplicate"],
                                      "genomes/hmq.bins.drep.out"),
            filtering_genome_min_length = config["params"]["dereplicate"]["drep"]["filtering_genome_min_length"],
            filtering_completeness = config["params"]["dereplicate"]["drep"]["filtering_completeness"],
            filtering_contamination = config["params"]["dereplicate"]["drep"]["filtering_contamination"],
            genome_comparison_algorithm = config["params"]["dereplicate"]["drep"]["genome_comparison_algorithm"],
            clustering_primary_ANI = config["params"]["dereplicate"]["drep"]["clustering_primary_ANI"],
            clustering_secondary_ANI = config["params"]["dereplicate"]["drep"]["clustering_secondary_ANI"],
            cov_thresh = config["params"]["dereplicate"]["drep"]["cov_thresh"],
            coverage_method = config["params"]["dereplicate"]["drep"]["coverage_method"],
            cluster_algorithm = config["params"]["dereplicate"]["drep"]["cluster_algorithm"],
            external_params = config["params"]["dereplicate"]["drep"]["external_params"]
        threads:
            config["params"]["dereplicate"]["drep"]["threads"]
        shell:
            '''
            dRep dereplicate \
            --processors {threads} \
            --genomes {input.bins_hmq} \
            --genomeInfo {input.genomes_info} \
            --length {params.filtering_genome_min_length} \
            --completeness {params.filtering_completeness} \
            --contamination {params.filtering_contamination} \
            --S_algorithm {params.genome_comparison_algorithm} \
            --P_ani {params.clustering_primary_ANI} \
            --S_ani {params.clustering_secondary_ANI} \
            --cov_thresh {params.cov_thresh} \
            --coverage_method {params.coverage_method} \
            --clusterAlg {params.cluster_algorithm} \
            {params.external_params} \
            {params.output_dir} \
            2> {log}

            touch {output}
            '''


    rule dereplicate_mags_drep_report:
        input:
            genomes_info = os.path.join(config["output"]["dereplicate"],
                                        "genomes_info/checkm_table_genomes_info.all.tsv"),
            drep_done = os.path.join(config["output"]["dereplicate"],
                                     "genomes/hmq.bins.drep.out/drep_done")
        output:
            rep_genomes_info = os.path.join(config["output"]["dereplicate"],
                         "report/checkm_table_genomes_info.derep.tsv")
        run:
            import pandas as pd
            from glob import glob

            genomes_info = pd.read_csv(input.genomes_info, sep="\t")
            
            rep_dir = os.path.join(os.path.dirname(input.drep_done), "dereplicated_genomes")
            rep_list = [os.path.basename(i) for i in sorted(glob(f'''{rep_dir}/*.fa'''))]

            rep_df = pd.DataFrame({"genome": rep_list})
            rep_df_info = rep_df.merge(genomes_info, how="left", on="genome")

            rep_df_info.to_csv(output.rep_genomes_info, sep="\t", index=False)


    rule dereplicate_mags_drep_all:
        input:
            os.path.join(config["output"]["dereplicate"],
                         "genomes_info/checkm_table_genomes_info.all.tsv"),
            os.path.join(config["output"]["dereplicate"],
                         "genomes_info/bins_hmq.all.tsv"),
            os.path.join(config["output"]["dereplicate"],
                         "genomes/hmq.bins.drep.out/drep_done"),
            os.path.join(config["output"]["dereplicate"],
                         "report/checkm_table_genomes_info.derep.tsv")
 
else:
    rule dereplicate_mags_drep_all:
        input:


rule dereplicate_mags_all:
    input:
        rules.dereplicate_mags_drep_all.input


localrules:
    dereplicate_mags_prepare,
    dereplicate_mags_drep_all,
    dereplicate_mags_drep_report,
    dereplicate_mags_all