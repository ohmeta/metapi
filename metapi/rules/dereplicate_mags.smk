if config["params"]["checkm"]["do"]:
    rule dereplicate_mags_prepare:
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
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.tsv"),
            genomes_info_simple = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.simple.csv"),
            mags_hmq = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/MAGs_hmq.{assembler}.all.tsv")
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


    localrules:
        dereplicate_mags_prepare


if config["params"]["dereplicate"]["drep"]["do"] and config["params"]["checkm"]["do"]:
    rule dereplicate_mags_drep:
        input:
            genomes_info = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.tsv"),
            mags_hmq = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/MAGs_hmq.{assembler}.all.tsv")
        output:
            os.path.join(config["output"]["dereplicate"], "genomes/bacteriome/MAGs_hmq.{assembler}.drep.out/drep_done")
        log:
            os.path.join(config["output"]["dereplicate"], "logs/drep/MAGs_hmq.{assembler}.drep.log")
        benchmark:
            os.path.join(config["output"]["dereplicate"], "benchmark/drep/MAGs_hmq.{assembler}.drep.benchmark.txt")
        conda:
            config["envs"]["drep"]
        params:
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
            outdir=$(dirname {output})
            rm -rf $outdir

            dRep dereplicate \
            --processors {threads} \
            --genomes {input.mags_hmq} \
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
            $outdir \
            2> {log}

            touch {output}
            '''


    rule dereplicate_mags_drep_report:
        input:
            genomes_info = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.tsv"),
            drep_done = os.path.join(
                config["output"]["dereplicate"],
                "genomes/bacteriome/MAGs_hmq.{assembler}.drep.out/drep_done")
        output:
            rep_genomes_info = os.path.join(
                config["output"]["dereplicate"],
                "report/bacteriome/checkm_table_genomes_info.{assembler}.drep.tsv.gz")
        run:
            import pandas as pd
            from glob import glob

            genomes_info = pd.read_csv(input.genomes_info, sep="\t")
            
            rep_dir = os.path.join(os.path.dirname(input.drep_done), "dereplicated_genomes")
            rep_list = [os.path.basename(i) for i in sorted(glob(f'''{rep_dir}/*.fa.gz'''))]

            rep_df = pd.DataFrame({"genome": rep_list})
            rep_df_info = rep_df.merge(genomes_info, how="left", on="genome")

            rep_df_info.to_csv(output.rep_genomes_info, sep="\t", index=False)


    localrules:
        dereplicate_mags_drep_report


    rule dereplicate_mags_drep_all:
        input:
            expand([
                os.path.join(config["output"]["dereplicate"],
                             "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.tsv"),
                os.path.join(config["output"]["dereplicate"],
                             "genomes_info/bacteriome/MAGs_hmq.{assembler}.all.tsv"),
                os.path.join(config["output"]["dereplicate"],
                             "genomes/bacteriome/MAGs_hmq.{assembler}.drep.out/drep_done"),
                os.path.join(config["output"]["dereplicate"],
                             "report/bacteriome/checkm_table_genomes_info.{assembler}.drep.tsv.gz")],
                assembler=ASSEMBLERS)
 
else:
    rule dereplicate_mags_drep_all:
        input:


if config["params"]["dereplicate"]["galah"]["do"] and config["params"]["checkm"]["do"]:
    rule dereplicate_mags_galah:
        input:
            genomes_info = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.simple.csv"),
            mags_hmq = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/MAGs_hmq.{assembler}.all.tsv")
        output:
            os.path.join(config["output"]["dereplicate"], "genomes/bacteriome/MAGs_hmq.{assembler}.galah.out/galah_done")
        log:
            os.path.join(config["output"]["dereplicate"], "logs/galah/MAGs_hmq.{assembler}.galah.log")
        benchmark:
            os.path.join(config["output"]["dereplicate"], "benchmark/galah/MAGs_hmq.{assembler}.galah.benchmark.txt")
        threads:
            config["params"]["dereplicate"]["galah"]["threads"]
        conda:
            config["envs"]["galah"]
        params:
            min_completeness = config["params"]["dereplicate"]["galah"]["min_completeness"],
            max_contamination = config["params"]["dereplicate"]["galah"]["max_contamination"],
            ani = config["params"]["dereplicate"]["galah"]["ani"],
            min_aligned_fraction = config["params"]["dereplicate"]["galah"]["min_aligned_fraction"],
            fragment_length = config["params"]["dereplicate"]["galah"]["fragment_length"],
            quality_formula = config["params"]["dereplicate"]["galah"]["quality_formula"],
            precluster_ani = config["params"]["dereplicate"]["galah"]["precluster_ani"],
            precluster_method = config["params"]["dereplicate"]["galah"]["precluster_method"],
            output_cluster_definition = os.path.join(
                config["output"]["dereplicate"],
                "genomes/bacteriome/MAGs_hmq.{assembler}.galah.out",
                "cluster_definition.tsv"),
            output_representative_list = os.path.join(
                config["output"]["dereplicate"],
                "genomes/bacteriome/MAGs_hmq.{assembler}.galah.out",
                "cluster_representative.list"),
            genome_dir = os.path.join(
                config["output"]["dereplicate"],
                "genomes/bacteriome/MAGs_hmq.{assembler}.galah.out",
                "cluster_representative.list"),
 
        shell:
            '''
            outdir=$(dirname {output})
            rm -rf $outdir
            mkdir -p $outdir

            galah cluster \
            --genome-fasta-list {input.mags_hmq} \
            --genome-info {input.genomes_info} \
            --genome-fasta-extension gz \
            --min-completeness {params.min_completeness} \
            --max-contamination {params.max_contamination} \
            --ani {params.ani} \
            --min-aligned-fraction {params.min_aligned_fraction} \
            --fragment-length {params.fragment_length} \
            --quality-formula {params.quality_formula} \
            --precluster-ani {params.precluster_ani} \
            --precluster-method {params.precluster_method} \
            --output-cluster-definition {params.output_cluster_definition} \
            --output-representative-list {params.output_representative_list} \
            --output-representative-fasta-directory $outdir/dereplicated_genomes \
            --threads {threads} \
            > {log} 2>&1

            touch {output}
            '''


    rule dereplicate_mags_galah_report:
        input:
            genomes_info = os.path.join(
                config["output"]["dereplicate"],
                "genomes_info/bacteriome/checkm_table_genomes_info.{assembler}.all.tsv"),
            galah_done = os.path.join(
                config["output"]["dereplicate"],
                "genomes/bacteriome/MAGs_hmq.{assembler}.galah.out/galah_done")
        output:
            rep_genomes_info = os.path.join(
                config["output"]["dereplicate"],
                "report/bacteriome/checkm_table_genomes_info.{assembler}.galah.tsv.gz")
        run:
            import pandas as pd
            from glob import glob

            genomes_info = pd.read_csv(input.genomes_info, sep="\t")
            
            rep_dir = os.path.join(os.path.dirname(input.galah_done), "dereplicated_genomes")
            rep_list = [os.path.basename(i) for i in sorted(glob(f'''{rep_dir}/*.fa.gz'''))]

            rep_df = pd.DataFrame({"genome": rep_list})
            rep_df_info = rep_df.merge(genomes_info, how="left", on="genome")

            rep_df_info.to_csv(output.rep_genomes_info, sep="\t", index=False)


    localrules:
        dereplicate_mags_galah_report


    rule dereplicate_mags_galah_all:
        input:
            expand([
                os.path.join(config["output"]["dereplicate"],
                             "genomes/bacteriome/MAGs_hmq.{assembler}.galah.out/galah_done"),
                os.path.join(config["output"]["dereplicate"],
                             "report/bacteriome/checkm_table_genomes_info.{assembler}.galah.tsv.gz")],
                assembler=ASSEMBLERS)

else:
    rule dereplicate_mags_galah_all:
        input:
 

rule dereplicate_mags_all:
    input:
        rules.dereplicate_mags_drep_all.input,
        rules.dereplicate_mags_galah_all.input


localrules:
    dereplicate_mags_drep_all,
    dereplicate_mags_galah_all,
    dereplicate_mags_all