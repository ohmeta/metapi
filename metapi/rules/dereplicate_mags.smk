if config["params"]["dereplicate"]["drep"]["do"]:
    rule dereplicate_mags_drep_prepare:
        input:
            bins_report = os.path.join(
                config["output"]["binning"],
                "report/assembly_stats_{assembler}_{binner_checkm}.tsv"),
            checkm_table = os.path.join(
                config["output"]["checkm"],
                "report/checkm_table_{assembler}_{binner_checkm}.tsv")
        output:
            genome_info = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/genomes_info_{assembler}_{binner_checkm}_checkm.csv")
        run:
            import pandas as pd

            bins_report = pd.read_csv(input.bins_report, sep='\t', header=[0, 1])\
                            .rename(columns={
                                "Unnamed: 0_level_1": "assembly_group",
                                "Unnamed: 1_level_1": "bin_id",
                                "Unnamed: 2_level_1": "bin_file",
                                "Unnamed: 3_level_1": "assembler",
                                "Unnamed: 4_level_1": "binner",
                            }, level=1)
            bins_report = bins_report[[
                ("assembly_group", "assembly_group"),
                ("bin_id", "bin_id"),
                ("bin_file", "bin_file"),
                ("assembler", "assembler"),
                ("binner", "binner"),
                ("length", "sum"),
                ("length", "N50")]]
            bins_report.columns = ["assembly_group", "bin_id", "bin_file", "assembler", "binner", "length", "N50"]
            bins_report = bins_report.sort_values(by=["bin_id"]).set_index("bin_id")

            checkm_table = pd.read_csv(input.checkm_table, sep='\t')\
                             .sort_values(by=["bin_id"]).set_index("bin_id")

            genome_info = pd.concat([bins_report, checkm_table], axis=1)\
                            .reset_index()\
                            .rename(columns={"index": "bin_id"})
            genome_info = genome_info.rename(columns={"bin_file": "genome"})
            genome_info.to_csv(output.genome_info, index=False)


    rule dereplicate_mags_drep:
        input:
            bins_hmq = os.path.join(config["output"]["checkm"],
                                    "report/{assembler}_{binner_checkm}_bins_hmq.tsv"),
            genome_info = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/genomes_info_{assembler}_{binner_checkm}_checkm.csv")
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
            {params.output_dir} \
            --processors {threads} \
            --length {params.filtering_genome_min_length} \
            --completeness {params.filtering_completeness} \
            --contamination {params.filtering_contamination} \
            --S_algorithm {params.genome_comparison_algorithm} \
            --P_ani {params.clustering_primary_ANI} \
            --S_ani {params.clustering_secondary_ANI} \
            --genomes {input.bins_hmq} \
            --genomeInfo {input.genome_info} \
            2> {log}

            touch {output}
            '''


    rule dereplicate_mags_drep_all:
        input:
            expand([
                os.path.join(
                    config["output"]["dereplicate"],
                    "genomes_info/genomes_info_{assembler}_{binner_checkm}_checkm.csv"),
                os.path.join(
                    config["output"]["dereplicate"],
                    "genomes/hmq.bins.{assembler}.{binner_checkm}.drep.out/drep_done")],
                   assembler=ASSEMBLERS,
                   binner_checkm=BINNERS_CHECKM)

else:
    rule dereplicate_mags_drep_all:
        input:


rule dereplicate_mags_all:
    input:
        rules.dereplicate_mags_drep_all.input