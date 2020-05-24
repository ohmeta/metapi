if config["params"]["dereplicate"]["drep"]["do"]:
    rule dereplicate_drep:
        input:
            bins_hmq = os.path.join(
                config["output"]["checkm"],
                "bins_hmq/{assembler}.{binner}.links"),
            genome_info = os.path.join(
                config["output"]["checkm"],
                "report/{assembler}_{binner}_checkm_table.tsv")
        output:
            directory(os.path.join(
                config["output"]["dereplicate"],
                "hmq.bins.{assembler}.{binner}.drep.out"))
        log:
            os.path.join(config["output"]["dereplicate"],
                         "logs/hmq.bins.{assembler}.{binner}.drep.log")
        params:
            bin_suffix = "fa",
            filtering_genome_min_length = \
                config["params"]["dereplicate"]["drep"]["filtering_genome_min_length"],
            filtering_completeness = \
                config["params"]["dereplicate"]["drep"]["filtering_completeness"],
            filtering_contamination = \
                config["params"]["dereplicate"]["drep"]["filtering_contamination"],
            genome_comparison_algorithm = \
                config["params"]["dereplicate"]["drep"]["genome_comparison_algorithm"],
            clustering_primary_ANI = \
                config["params"]["dereplicate"]["drep"]["clustering_primary_ANI"],
            clustering_secondary_ANI = \
                config["params"]["dereplicate"]["drep"]["clustering_secondary_ANI"]
        threads:
            config["params"]["dereplicate"]["drep"]["threads"]
        shell:
            '''
            dRep dereplicate \
            {output} \
            --processors {threads} \
            --length {params.filtering_genome_min_length} \
            --completeness {params.filtering_completeness} \
            --contamination {params.filtering_contamination} \
            --S_algorithm {params.genome_comparison_algorithm} \
            --P_ani {params.clustering_primary_ANI} \
            --S_ani {params.clustering_secondary_ANI} \
            --genomes {input.bins_hmq}/*.{params.bin_suffix} \
            --genomeInfo {input.genome_info} \
            2> {log}
            '''


    rule dereplicate_drep_all:
        input:
            expand(
                os.path.join(config["output"]["dereplicate"],
                             "hmq.bins.{assembler}.{binner}.drep.out"),
                assembler=ASSEMBLERS,
                binner=BINNERS)

else:
    rule dereplicate_drep_all:
        input:


rule dereplicate_all:
    input:
        rules.dereplicate_drep_all.input
