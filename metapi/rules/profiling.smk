rule metaphlan2_profiling:
    input:
        reads = clean_reads
    output:
        bt2_out = os.path.join(config["results"]["profiling"]["metaphlan2"]["bowtie2_out"], "{sample}.bowtie2.bz2"),
        profile = os.path.join(config["results"]["profiling"]["metaphlan2"]["profile"], "{sample}.metaphlan2.profile")
    params:
        bowtie2db = config["params"]["profiling"]["metaphlan2"]["bowtie2db"],
        index = config["params"]["profiling"]["metaphlan2"]["index"],
        input_type = config["params"]["profiling"]["metaphlan2"]["input_type"],
        taxonomic_level = config["params"]["profiling"]["metaphlan2"]["taxonomic_level"],
        analysis_type = config["params"]["profiling"]["metaphlan2"]["analysis_type"],
        min_cu_len = config["params"]["profiling"]["metaphlan2"]["min_cu_len"],
        read_min_len = config["params"]["profiling"]["metaphlan2"]["read_min_len"],
        no_unknown_estimation = "--no_unknown_estimation" if config["params"]["profiling"]["metaphlan2"]["no_unknown_estimation"] else "",
        sample_id = "{sample}"
    log:
        os.path.join(config["logs"]["profiling"]["metaphlan2"], "{sample}_metaphlan2.log")
    threads:
        config["params"]["profiling"]["metaphlan2"]["threads"]
    shell:
        '''
        metaphlan2.py \
        {input.reads[0]},{input.reads[1]} \
        --nproc {threads} \
        --min_cu_len {params.min_cu_len} \
        --bowtie2db {params.bowtie2db} \
        --index {params.index} \
        --input_type {params.input_type} \
        --tax_lev {params.taxonomic_level} \
        -t {params.analysis_type} \
        {params.no_unknown_estimation} \
        --sample_id {params.sample_id} \
        --bowtie2out {output.bt2_out} \
        > {output.profile} \
        2> {log}
        '''

rule metaphlan2_merge:
    input:
        expand("{profile}/{sample}.metaphlan2.profile",
               profile=config["results"]["profiling"]["metaphlan2"]["profile"],
               sample=_samples.index.unique())
    output:
        os.path.join(config["results"]["profiling"]["metaphlan2"]["base_dir"], "metaphlan2.merged.profile")
    log:
        os.path.join(config["logs"]["profiling"]["metaphlan2"], "metaphlan2.merged.log")
    shell:
        '''
        merge_metaphlan_tables.py {input} > {output} 2> {log}
        sed -i 's/.metaphlan2//g' {output}
        '''


rule jgi_profiling:
    input:
        reads = clean_reads,
        index_database = expand("{prefix}.{suffix}",
                                prefix=config["params"]["profiling"]["jgi"]["index_prefix"],
                                suffix=["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l", "rev.1.bt2l", "rev.2.bt2l"])
    output:
        depth = os.path.join(config["results"]["profiling"]["jgi"]["depth"], "{sample}.jgi.depth.gz")
    log:
        os.path.join(config["logs"]["profiling"]["jgi"], "{sample}.profiling.jgi.log")
    threads:
        config["params"]["profiling"]["threads"]
    params:
        index_prefix = config["params"]["profiling"]["jgi"]["index_prefix"],
        memory_limit = config["params"]["profiling"]["jgi"]["memory_limit"],
        fragment = config["params"]["profiling"]["jgi"]["fragment"],
        temp_prefix = os.path.join(config["results"]["profiling"]["jgi"]["depth"], "temp.{sample}"),
        depth = os.path.join(config["results"]["profiling"]["jgi"]["depth"], "{sample}.jgi.depth"),
    run:
        if IS_PE:
            shell('''bowtie2 -x {params.index_prefix} -1 {input.reads[0]} -2 {input.reads[1]} \
                  --end-to-end --very-sensitive --phred33 --threads {threads} \
                  --seed 0 --time -k 2 --no-unal --no-discordant \
                  -X {params.fragment} 2> {log} | \
                  samtools sort -@{threads} -m {params.memory_limit} -T {params.temp_prefix} -O BAM - | \
                  jgi_summarize_bam_contig_depths --outputDepth {params.depth} -''')
        else:
            shell('''bowtie2 -x {params.index_prefix} -U {input.reads[0]} \
                  --end-to-end --very-sensitive --phred33 --threads {threads} \
                  --seed 0 --time -k 2 --no-unal --no-discordant \
                  -X {params.fragment} 2> {log} | \
                  samtools sort -@{threads} -m {params.memory_limit} -T {params.temp_prefix} -O BAM - | \
                  jgi_summarize_bam_contig_depths --outputDepth {params.depth} -''')

        shell('''pigz {params.depth}''')
        shell('''rm -rf {params.temp_prefix}*.bam''')
        shell('''echo "profiling done" >> {log}''')


rule jgi_profile_merge:
    input:
        abun_files = expand(os.path.join(config["results"]["profiling"]["jgi"]["depth"],
                                         "{sample}.jgi.depth.gz"),
                            sample=_samples.index.unique()),
        taxonomy = config["params"]["profiling"]["taxonomy"],
        index_metadata = config["params"]["profiling"]["index_metadata"],
    output:
        abundance_profile = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                         "abundance_profile.tsv"),
        depth_profile = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                     "depth_profile.tsv"),
        abundance_profile_k = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_superkingdom.tsv"),
        abundance_profile_p = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_phylum.tsv"),
        abundance_profile_o = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_order.tsv"),
        abundance_profile_c = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_class.tsv"),
        abundance_profile_f = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_family.tsv"),
        abundance_profile_g = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_genus.tsv"),
        abundance_profile_s = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_species.tsv"),
        abundance_profile_t = os.path.join(config["results"]["profiling"]["jgi"]["profile"],
                                           "abundance_profile_strain.tsv")
    threads:
        config["params"]["profiling"]["threads"]
    run:
        import pandas as pd
        taxonomy_df = pd.read_csv(input.taxonomy, sep='\t')
        merger.global_init(input.index_metadata)
        depth_df, abun_df = merger.get_all_abun_df(input.abun_files, threads, "jgi")
        depth_df.to_csv(output.depth_profile, sep='\t', index=False)
        abun_df.to_csv(output.abundance_profile, sep='\t', index=False)
        samples_list = sorted(abun_df.columns[1:].to_list())
        abun_tax_df = abun_df.merge(taxonomy_df)

        merger.get_profile(abun_tax_df, samples_list,  "lineages_superkingdom_new", output.abundance_profile_k)
        merger.get_profile(abun_tax_df, samples_list, "lineages_phylum_new", output.abundance_profile_p)
        merger.get_profile(abun_tax_df, samples_list, "lineages_order_new", output.abundance_profile_o)
        merger.get_profile(abun_tax_df, samples_list, "lineages_class_new", output.abundance_profile_c)
        merger.get_profile(abun_tax_df, samples_list, "lineages_family_new", output.abundance_profile_f)
        merger.get_profile(abun_tax_df, samples_list, "lineages_genus_new", output.abundance_profile_g)
        merger.get_profile(abun_tax_df, samples_list, "lineages_species_new", output.abundance_profile_s)
        merger.get_profile(abun_tax_df, samples_list, "lineages_strain_new", output.abundance_profile_t)


rule humann2_profiling:
    input:
        reads = clean_reads
    output:
        reads_merged = temp(os.path.join(config["results"]["profiling"]["humann2"],
                                         "{sample}.humann2_out/{sample}.merged.fq.gz")),
        genefamilies = os.path.join(config["results"]["profiling"]["humann2"],
                                    "{sample}.humann2_out/{sample}_genefamilies.tsv"),
        pathabundance = os.path.join(config["results"]["profiling"]["humann2"],
                                     "{sample}.humann2_out/{sample}_pathabundance.tsv"),
        pathcoverage = os.path.join(config["results"]["profiling"]["humann2"],
                                    "{sample}.humann2_out/{sample}_pathcoverage.tsv")
    params:
        base_name = "{sample}",
        remove = "--remove-temp-output" if config["params"]["profiling"]["humann2"]["remove-temp-output"] else "",
        out_dir = os.path.join(config["results"]["profiling"]["humann2"], "{sample}.humann2_out")
    threads:
        8
    log:
        os.path.join(config["logs"]["profiling"]["humann2"], "{sample}_humann2.log")
    run:
        reads = input.reads
        if IS_PE:
            reads_str = " ".join(input.reads)
            shell('''cat %s > %s''' % (reads_str, output.reads_merged))
            reads = output.reads_merged

        shell('''humann2 --input %s --output %s --output-basename %s --threads %d --o-log %s %s''' % \
              (reads, params.out_dir, params.base_name, threads, log, params.remove))


rule humann2_postprocess:
    input:
        genefamilies = os.path.join(config["results"]["profiling"]["humann2"],
                                    "{sample}.humann2_out/{sample}_genefamilies.tsv"),
        pathabundance = os.path.join(config["results"]["profiling"]["humann2"],
                                     "{sample}.humann2_out/{sample}_pathabundance.tsv"),
        pathcoverage = os.path.join(config["results"]["profiling"]["humann2"],
                                    "{sample}.humann2_out/{sample}_pathcoverage.tsv")
    output:
        genefamilies = expand("{humann2}/{{sample}}.humann2_out/{{sample}}_genefamilies.{norm}.tsv",
                              humann2=config["results"]["profiling"]["humann2"],
                              norm=config["params"]["profiling"]["humann2"]["normalize_method"]),
        pathabundance = expand("{humann2}/{{sample}}.humann2_out/{{sample}}_pathabundance.{norm}.tsv",
                               humann2=config["results"]["profiling"]["humann2"],
                               norm=config["params"]["profiling"]["humann2"]["normalize_method"]),
        pathcoverage = expand("{humann2}/{{sample}}.humann2_out/{{sample}}_pathcoverage.{norm}.tsv",
                              humann2=config["results"]["profiling"]["humann2"],
                              norm=config["params"]["profiling"]["humann2"]["normalize_method"])
        groupprofile = expand("{humann2}/{{sample}}.humann2_out/{{sample}}_group-{group}-profile.tsv",
                              humann2=config["results"]["profiling"]["humann2"],
                              group=config["params"]["profiling"]["humann2"]["map_database"])

    params:
        normalize_method = config["params"]["profiling"]["humann2"]["normalize_method"],
        regroup_method = config["params"]["profiling"]["humann2"]["regroup_method"],
        map_database =  config["params"]["profiling"]["humann2"]["map_database"]
    run:
        shell("humann2_renorm_table --input {input.genefamilies} --update-snames \
              --output {output.genefamilies} --units {params.normalize_method}")
        shell("humann2_renorm_table --input {input.pathabundance} --update-snames \
              --output {output.pathabundance} --units {params.normalize_method}")
        shell("humann2_renorm_table --input {input.pathcoverage} --update-snames \
              --output {output.pathcoverage} --units {params.normalize_method}")
        i = -1
        for db in params.map_database:
            i += 1
            shell("humann2_regroup_table --input {input.genefamilies} --groups %s \
                  --function {params.regroup_method} --output {output.groupprofile[i]}" % db)
