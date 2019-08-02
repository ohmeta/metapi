rule metaphlan2_profilling:
    input:
        reads = clean_reads
    output:
        bt2_out = os.path.join(config["results"]["profilling"]["metaphlan2"]["bowtie2_out"], "{sample}.bowtie2.bz2"),
        profile = os.path.join(config["results"]["profilling"]["metaphlan2"]["profile"], "{sample}.metaphlan2.profile")
    params:
        metaphlan2_env = config["params"]["profilling"]["metaphlan2"]["env"],
        metaphlan2_script = config["params"]["profilling"]["metaphlan2"]["script"],
        bowtie2_exe = config["params"]["profilling"]["metaphlan2"]["bowtie2_exe"],
        input_type = config["params"]["profilling"]["metaphlan2"]["input_type"],
        mpa_pkl = config["params"]["profilling"]["metaphlan2"]["mpa_pkl"],
        bowtie2db = config["params"]["profilling"]["metaphlan2"]["bowtie2db"],
        min_cu_len = config["params"]["profilling"]["metaphlan2"]["min_cu_len"],
        sample_id = "{sample}"
    log:
        os.path.join(config["logs"]["profilling"]["metaphlan2"], "{sample}_metaphlan2.log")
    threads:
        config["params"]["profilling"]["metaphlan2"]["threads"]
    shell:
        '''
        #set +u; source activate {params.metaphlan2_env}; set -u;
        /ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv2/bin/python \
        {params.metaphlan2_script} \
        {input.reads[0]},{input.reads[1]} \
        --bowtie2_exe {params.bowtie2_exe} \
        --mpa_pkl {params.mpa_pkl} \
        --bowtie2db {params.bowtie2db} \
        --min_cu_len {params.min_cu_len} \
        --bowtie2out {output.bt2_out} \
        --nproc {threads} \
        --input_type {params.input_type} \
        --sample_id {params.sample_id} \
        > {output.profile} \
        2> {log}
        '''

rule metaphlan2_merge:
    input:
        expand("{profile}/{sample}.metaphlan2.profile",
               profile=config["results"]["profilling"]["metaphlan2"]["profile"],
               sample=_samples.index.unique())
    output:
        os.path.join(config["results"]["profilling"]["metaphlan2"]["base_dir"], "metaphlan2.merged.profile")
    log:
        os.path.join(config["logs"]["profilling"]["metaphlan2"], "metaphlan2.merged.log")
    params:
        metaphlan2_env = config["params"]["profilling"]["metaphlan2"]["env"]
    shell:
        '''
        #!/bin/bash
        #set +u; source activate {params.metaphlan2_env}; set -u;
        /ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv2/bin/python /ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv2/bin/merge_metaphlan_tables.py {input} > {output} 2> {log}
        sed -i 's/.metaphlan2//g' {output}
        '''

rule mwas_profilling:
    input:
        reads = clean_reads,
        index_metadata = config["params"]["profilling"]["index_metadata"],
        index_database = expand("{prefix}.{suffix}",
                                prefix=config["params"]["profilling"]["comg"]["index_prefix"],
                                suffix=["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l", "rev.1.bt2l", "rev.2.bt2l"])
    output:
        abundance= os.path.join(config["results"]["profilling"]["comg"]["abundance"], "{sample}.comg.abundance.gz"),
        depth = os.path.join(config["results"]["profilling"]["metabat2"]["depth"], "{sample}.metabat2.depth.gz")
    log:
        os.path.join(config["logs"]["profilling"]["comg"], "{sample}.profilling.comg.log")
    threads:
        config["params"]["profilling"]["threads"]
    params:
        sample = "{sample}",
        index_prefix = config["params"]["profilling"]["comg"]["index_prefix"],
        sam_format = config["params"]["profilling"]["sam_format"],
        identity = config["params"]["profilling"]["comg"]["identity"],
        fragment = config["params"]["profilling"]["comg"]["fragment"],
        output_type = config["params"]["profilling"]["comg"]["output_type"],
        insert_size = config["params"]["profilling"]["comg"]["insert_size"],
        abundance_script = config["params"]["profilling"]["comg"]["abundance_script"],
        abundance_outdir = config["results"]["profilling"]["comg"]["abundance"],
        samtools_sort_prefix = os.path.join(config["results"]["profilling"]["comg"]["abundance"], "{sample}.temp"),
        depth = os.path.join(config["results"]["profilling"]["metabat2"]["depth"], "{sample}.metabat2.depth"),
    run:
        if IS_PE:
            shell('''bowtie2 -x {params.index_prefix} -1 {input.reads[0]} -2 {input.reads[1]} \
                  --end-to-end --very-sensitive --phred33 --threads {threads} \
                  --seed 0 --time -k 2 --no-unal --no-discordant \
                  -X {params.fragment} 2> {log} | \
                  samtools sort -@{threads} -T {params.samtools_sort_prefix} -O {params.sam_format} - | \
                  tee >(python {params.abundance_script} \
                  --sample-name {params.sample} \
                  --sam-file - \
                  --sam-format {params.sam_format} \
                  --insert-size {params.insert_size} \
                  --marker-matrix {input.index_metadata} \
                  --outdir {params.abundance_outdir} \
                  --identity {params.identity} \
                  --output-type {params.output_type}) | \
                  jgi_summarize_bam_contig_depths --outputDepth {params.depth} -''')
        else:
            shell('''bowtie2 -x {params.index_prefix} -U {input.reads[0]} \
                  --end-to-end --very-sensitive --phred33 --threads {threads} \
                  --seed 0 --time -k 2 --no-unal --no-discordant \
                  -X {params.fragment} 2> {log} | \
                  samtools sort -@{threads} -T {params.samtools_sort_prefix} -O {params.sam_format} - | \
                  tee >(python {params.abundance_script} \
                  --sample-name {params.sample} \
                  --sam-file - \
                  --sam-format {params.sam_format} \
                  --insert-size {params.insert_size} \
                  --marker-matrix {input.index_metadata} \
                  --outdir {params.abundance_outdir} \
                  --identity {params.identity} \
                  --output-type {params.output_type}) | \
                  jgi_summarize_bam_contig_depths --outputDepth {params.depth} -''')

        shell('''pigz {params.depth}''')
        shell('''pigz -c {params.abundance_outdir}/{params.sample}.abundance > {output.abundance}''')
        shell('''rm -rf {params.abundance_outdir}/{params.sample}.abundance''')
        shell('''rm -rf {params.samtools_sort_prefix}.*.bam''')
        shell('''echo "profilling done" >> {log}''')


rule mwas_profile_merge_hsx:
    input:
        abun_files = expand(os.path.join(config["results"]["profilling"]["comg"]["abundance"], "{sample}.comg.abundance.gz"),
                            sample=_samples.index.unique()),
        taxonomy = config["params"]["profilling"]["taxonomy"]
    output:
        abundance_profile = os.path.join(config["results"]["profilling"]["comg"]["profile"], "abundance_profile.tsv"),
        count_profile = os.path.join(config["results"]["profilling"]["comg"]["profile"], "count_profile.tsv"),
        abundance_profile_k = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_superkingdom.tsv"),
        abundance_profile_p = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_phylum.tsv"),
        abundance_profile_o = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_order.tsv"),
        abundance_profile_c = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_class.tsv"),
        abundance_profile_f = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_family.tsv"),
        abundance_profile_g = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_genus.tsv"),
        abundance_profile_s = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_species.tsv"),
        abundance_profile_t = os.path.join(config["results"]["profilling"]["comg"]["profile"],
                                           "abundance_profile_strain.tsv")
    threads:
        config["params"]["profilling"]["threads"]
    run:
        import pandas as pd

        taxonomy_df = pd.read_csv(input.taxonomy, sep='\t')

        count_df, abun_df = merger.get_all_abun_df(input.abun_files, threads, "hsx")

        count_df.to_csv(output.count_profile, sep='\t', index=False)
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

rule mwas_profile_merge_jgi:
    input:
        abun_files = expand(os.path.join(config["results"]["profilling"]["metabat2"]["depth"], "{sample}.metabat2.depth.gz"),
                            sample=_samples.index.unique()),
        taxonomy = config["params"]["profilling"]["taxonomy"],
        index_metadata = config["params"]["profilling"]["index_metadata"],
    output:
        abundance_profile = os.path.join(config["results"]["profilling"]["metabat2"]["profile"], "abundance_profile.tsv"),
        depth_profile = os.path.join(config["results"]["profilling"]["metabat2"]["profile"], "depth_profile.tsv"),
        abundance_profile_k = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_superkingdom.tsv"),
        abundance_profile_p = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_phylum.tsv"),
        abundance_profile_o = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_order.tsv"),
        abundance_profile_c = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_class.tsv"),
        abundance_profile_f = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_family.tsv"),
        abundance_profile_g = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_genus.tsv"),
        abundance_profile_s = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_species.tsv"),
        abundance_profile_t = os.path.join(config["results"]["profilling"]["metabat2"]["profile"],
                                           "abundance_profile_strain.tsv")
    threads:
        config["params"]["profilling"]["threads"]
    run:
        import pandas as pd

        taxonomy_df = pd.read_csv(input.taxonomy, sep='\t')

        global INDEX_METADATA
        INDEX_METADATA = pd.read_csv(input.index_metadata, sep='\t')

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
