if config["params"]["dereplicate"]["drep"]["do"]:
    rule dereplicate_mags_drep:
        input:
            genomes_info = os.path.join(config["output"]["checkm"],
                                        "report/checkm_table_{assembler}_{binner_checkm}.tsv"),
            bins_hmq = os.path.join(config["output"]["checkm"],
                                    "report/{assembler}_{binner_checkm}_bins_hmq.tsv")
        output:
            os.path.join(config["output"]["dereplicate"],
                         "genomes/hmq.bins.{assembler}.{binner_checkm}.drep.out/drep_done")
        log:
            os.path.join(config["output"]["dereplicate"],
                         "logs/hmq.bins.{assembler}.{binner_checkm}.drep.log")
        benchmark:
            os.path.join(config["output"]["dereplicate"],
                         "benchmark/{assembler}.{binner_checkm}.drep.benchmark.txt")
        conda:
            config["envs"]["drep"]
        params:
            output_dir = os.path.join(config["output"]["dereplicate"],
                                      "genomes/hmq.bins.{assembler}.{binner_checkm}.drep.out"),
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


    rule dereplicate_mags_drep_all:
        input:
            expand(
                os.path.join(
                    config["output"]["dereplicate"],
                    "genomes/hmq.bins.{assembler}.{binner_checkm}.drep.out/drep_done"),
                   assembler=ASSEMBLERS,
                   binner_checkm=BINNERS_CHECKM)

else:
    rule dereplicate_mags_drep_all:
        input:


rule dereplicate_mags_all:
    input:
        rules.dereplicate_mags_drep_all.input