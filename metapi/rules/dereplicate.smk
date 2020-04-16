rule drep:
    input:
        bins_hmq = os.path.join(config["results"]["checkm"]["base_dir"],
                                "bins.{assembler}.{binner}_out.hmq"),
        genome_info = os.path.join(config["results"]["checkm"]["base_dir"],
                                   "{assembler}.{binner}.checkm.out.tsv")
    output:
        out_dir = directory(os.path.join(config["results"]["dereplication"],
                                         "hmq.bins.{assembler}.{binner}.drep_out"))
    params:
        filtering_genome_min_length = config["params"]["dereplication"]["drep"]["filtering_genome_min_length"],
        filtering_completeness = config["params"]["dereplication"]["drep"]["filtering_completeness"],
        filtering_contamination = config["params"]["dereplication"]["drep"]["filtering_contamination"],
        genome_comparison_algorithm = config["params"]["dereplication"]["drep"]["genome_comparison_algorithm"],
        clustering_primary_ANI = config["params"]["dereplication"]["drep"]["clustering_primary_ANI"],
        clustering_secondary_ANI = config["params"]["dereplication"]["drep"]["clustering_secondary_ANI"]
    log:
        os.path.join(config["logs"]["dereplication"], "{assembler}.{binner}.drep.log")
    threads:
        config["params"]["dereplication"]["drep"]["threads"]
    shell:
        '''
        rm -rf {output.out_dir}

        dRep dereplicate {output.out_dir} \
        --processors {threads} \
        --length {params.filtering_genome_min_length} \
        --completeness {params.filtering_completeness} \
        --contamination {params.filtering_contamination} \
        --S_algorithm {params.genome_comparison_algorithm} \
        --P_ani {params.clustering_primary_ANI} \
        --S_ani {params.clustering_secondary_ANI} \
        --genomes {input.bins_hmq}/*.fa \
        --genomeInfo {input.genome_info} 2> {log}
        '''
