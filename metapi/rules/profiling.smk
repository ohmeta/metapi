if config["params"]["profiling"]["metaphlan2"]["do"]:
    rule profiling_metaphlan2:
        input:
            assembly_input
        output:
           protected(os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"))
        log:
            os.path.join(
                config["output"]["profiling"],
                         "logs/{sample}.metaphlan2.log")
        params:
            sample_id = "{sample}",
            input_type = config["params"]["profiling"]["metaphlan2"]["input_type"],
            read_min_len = config["params"]["profiling"]["metaphlan2"]["read_min_len"],
            bowtie2db = config["params"]["profiling"]["metaphlan2"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan2"]["index"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan2"]["bowtie2_presets"],
            min_cu_len = config["params"]["profiling"]["metaphlan2"]["min_cu_len"],
            taxonomic_level = config["params"]["profiling"]["metaphlan2"]["taxonomic_level"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan2"]["avoid_disqm"] \
                   else "",
            stat_q = config["params"]["profiling"]["metaphlan2"]["stat_q"],
            stat = config["params"]["profiling"]["metaphlan2"]["stat"],
            analysis_type = config["params"]["profiling"]["metaphlan2"]["analysis_type"],
            no_unknown_estimation = "--no_unknown_estimation" \
                if config["params"]["profiling"]["metaphlan2"]["no_unknown_estimation"] \
                   else "",
            no_map = config["params"]["profiling"]["metaphlan2"]["no_map"],
            bowtie2_out = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.bowtie2.bz2"),
            biom = config["params"]["profiling"]["metaphlan2"]["biom"],
            biom_out = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.biom")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            shell(
                '''
                metaphlan2.py \
                %s \
                --input_type {params.input_type} \
                --read_min_len {params.read_min_len} \
                --nproc {threads} \
                %s \
                --bowtie2db {params.bowtie2db} \
                --index {params.index} \
                --bt2_ps {params.bowtie2_presets} \
                --min_cu_len {params.min_cu_len} \
                --tax_lev {params.taxonomic_level} \
                {params.avoid_disqm} \
                --stat_q {params.stat_q} \
                --stat {params.stat} \
                -t {params.analysis_type} \
                {params.no_unknown_estimation} \
                --output_file {output} \
                --sample_id {params.sample_id} \
                %s \
                2> {log}
                ''' % (
                    ",".join(input),

                    "--no_map" \
                    if params.no_map \
                    else "--bowtie2out %s" % params.bowtie2_out,

                    "--biom %s" % params.biom_out \
                    if params.biom \
                    else ""))


    rule profiling_metaphlan2_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"),
                   sample=SAMPLES.index.unique())
        output:
            os.path.join(
                config["output"]["profiling"],
                         "profile/metaphlan2.merged.abundance.profile.tsv")
        log:
            os.path.join(
                config["output"]["profiling"], "logs/metaphlan2.merged.log")
        shell:
            '''
            merge_metaphlan_tables.py \
            {input} \
            > {output} \
            2> {log}

            sed -i 's/.metaphlan2//g' {output}
            '''


    rule profiling_metaphlan2_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2.merged.abundance.profile.tsv")

else:
    rule profiling_metaphlan2_all:
        input:
           

if config["params"]["profiling"]["jgi"]["do"]:
    if not config["params"]["profiling"]["jgi"]["oneway"]:
        rule profiling_jgi:
            input:
                reads = assembly_input,
                index_database = expand(
                    "{prefix}.{suffix}",
                    prefix=config["params"]["profiling"]["jgi"]["index_prefix"],
                    suffix=["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l",
                            "rev.1.bt2l", "rev.2.bt2l"])
            output:
                protected(os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.jgi.coverage.gz"))
            log:
                os.path.join(
                    config["output"]["profiling"],
                    "logs/{sample}.jgi.log")
            threads:
                config["params"]["profiling"]["threads"]
            params:
                index_prefix = config["params"]["profiling"]["jgi"]["index_prefix"],
                memory_limit = config["params"]["profiling"]["jgi"]["memory_limit"],
                fragment = config["params"]["profiling"]["jgi"]["fragment"],
                temp_prefix = os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.temp"),
                coverage = os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.jgi.coverage")
            run:
                shell(
                    '''
                    bowtie2 \
                    -x {params.index_prefix} \
                    %s \
                    --end-to-end \
                    --very-sensitive \
                    --phred33 \
                    --threads {threads} \
                    --seed 0 \
                    --time \
                    -k 2 \
                    --no-unal \
                    --no-discordant \
                    -X {params.fragment} \
                    2> {log} | \
                    samtools sort \
                    -@{threads} \
                    -m {params.memory_limit} \
                    -T {params.temp_prefix} -O BAM - | \
                    jgi_summarize_bam_contig_depths \
                    --outputDepth {params.coverage} -
                    ''' % \
                    "-1 {input.reads[0]} -2 {input.reads[1]}" if IS_PE \
                    else "-U {input.reads[0]}")

                shell('''gzip {params.coverage}''')
                shell('''rm -rf {params.temp_prefix}*.bam''')
                shell('''echo "Profiling done" >> {log}''')

    else:
        rule profiling_jgi:
            input:
                reads = lambda wildcards: get_reads(wildcards, "raw"),
                index_database = expand(
                    "{prefix}.{suffix}",
                    prefix=config["params"]["profiling"]["jgi"]["index_prefix"],
                    suffix=["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l",
                            "rev.1.bt2l", "rev.2.bt2l"]),
                index_host = expand(
                    "{prefix}.{suffix}",
                    prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]) \
                    if RMHOST_DO else ""
            output:
                protected(os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.jgi.coverage.gz"))
            log:
                os.path.join(
                    config["output"]["profiling"],
                    "logs/{sample}.jgi.log")
            threads:
                config["params"]["profiling"]["threads"]
            params:
                adapter_trimming = '--disable_adapter_trimming' \
                    if config["params"]["trimming"]["fastp"]["disable_adapter_trimming"] else "",
                index_prefix = config["params"]["profiling"]["jgi"]["index_prefix"],
                host_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                memory_limit = config["params"]["profiling"]["jgi"]["memory_limit"],
                compression_level = config["params"]["profiling"]["jgi"]["compression_level"],
                fragment = config["params"]["profiling"]["jgi"]["fragment"],
                temp_dir = os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.temp"),
                coverage = os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.jgi.coverage")
            run:
                if TRIMMING_DO and RMHOST_DO:
                    shell('''date > {log}''')
                    shell(
                        '''
                        fastp %s {params.adapter_trimming} --thread {threads} \
                        --stdout --json /dev/null --html /dev/null 2>> {log} | \
                        bowtie2 --threads {threads} -x {params.host_prefix} %s - 2>> {log} | \
                        samtools fastq -@{threads} -N -f 12 -F 256 - |
                        bowtie2 \
                        -x {params.index_prefix} \
                        %s - \
                        --end-to-end \
                        --very-sensitive \
                        --phred33 \
                        --threads {threads} \
                        --seed 0 \
                        --time \
                        -k 2 \
                        --no-unal \
                        --no-discordant \
                        -X {params.fragment} \
                        2>> {log} | \
                        sambamba view -q --nthreads {threads} \
                        --compression-level {params.compression_level} \
                        --format bam \
                        --compression-level {params.compression_level} \
                        --sam-input /dev/stdin \
                        --output-filename /dev/stdout | \
                        sambamba sort -q --nthreads {threads} \
                        --memory-limit {params.memory_limit} \
                        --compression-level 0 \
                        --tmpdir {params.temp_dir} \
                        --out /dev/stdout /dev/stdin | \
                        jgi_summarize_bam_contig_depths \
                        --outputDepth {params.coverage} - \
                        2>> {log}
                        ''' % (
                            "--in1 {input.reads[0]} --in2 {input.reads[1]}" if IS_PE else "--in1 {input.reads[0]}",
                            "--interleaved" if IS_PE else "",
                            "--interleaved" if IS_PE else ""
                        )
                    )
                    shell('''pigz --processes {threads} {params.coverage}''')
                    shell('''rm -rf {params.temp_dir}''')
                    shell('''echo "jgi profiling done" >> {log}''')
                    shell('''date >> {log}''')

    rule profiling_jgi_merge:
        input:
            coverage = expand(os.path.join(
                config["output"]["profiling"],
                "profile/jgi/{sample}/{sample}.jgi.coverage.gz"),
                              sample=SAMPLES.index.unique()),
            taxonomy = config["params"]["profiling"]["jgi"]["taxonomy"],
            index_metadata = config["params"]["profiling"]["jgi"]["index_metadata"]
        output:
            abundance_profile = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.tsv"),
            coverage_profile = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.coverage.profile.tsv"),
            abundance_profile_k = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.superkingdom.tsv"),
            abundance_profile_p = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.phylum.tsv"),
            abundance_profile_o = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.order.tsv"),
            abundance_profile_c = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.class.tsv"),
            abundance_profile_f = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.family.tsv"),
            abundance_profile_g = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.genus.tsv"),
            abundance_profile_s = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.species.tsv"),
            abundance_profile_t = os.path.join(
                config["output"]["profiling"],
                "profile/jgi.merged.abundance.profile.strain.tsv")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            taxonomy_df = pd.read_csv(input.taxonomy, sep='\t')

            metapi.profiler_init(input.index_metadata)

            coverage_df, abundance_df = metapi.get_all_abun_df(
                input.coverage, threads, "jgi")

            samples_list = sorted(abundance_df.columns[1:].to_list())

            coverage_df.to_csv(output.coverage_profile, sep='\t', index=False)
            abundance_df.to_csv(output.abundance_profile, sep='\t', index=False)

            abundance_taxonomy_df = abundance_df.merge(taxonomy_df)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_superkingdom_new", output.abundance_profile_k)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_phylum_new", output.abundance_profile_p)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_order_new", output.abundance_profile_o)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_class_new", output.abundance_profile_c)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_family_new", output.abundance_profile_f)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_genus_new", output.abundance_profile_g)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_species_new", output.abundance_profile_s)

            metapi.get_profile(abundance_taxonomy_df, samples_list,
                               "lineages_strain_new", output.abundance_profile_t)


    rule profiling_jgi_all:
        input:
             expand([
                 os.path.join(
                     config["output"]["profiling"],
                     "profile/jgi.merged.{target}.profile.tsv"),
                 os.path.join(
                     config["output"]["profiling"],
                     "profile/jgi.merged.abundance.profile.{level}.tsv")],
                    target=["abundance", "coverage"],
                    level=["superkingdom",
                           "phylum",
                           "order",
                           "class",
                           "family",
                           "genus",
                           "species",
                           "strain"])

else:
    rule profiling_jgi_all:
        input:

           
if config["params"]["profiling"]["humann2"]["do"]:
    rule profiling_humann2:
        input:
            reads = assembly_input
        output:
            reads_merged = temp(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}.merged.fq.gz")),
            genefamilies = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_genefamilies.tsv")),
            pathabundance = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_pathabundance.tsv")),
            pathcoverage = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_pathcoverage.tsv"))
        params:
            base_name = "{sample}",
            remove = "--remove-temp-output" \
                if config["params"]["profiling"]["humann2"]["remove_temp_output"] \
                   else "",
            output_dir = os.path.join(
                config["output"]["profiling"],
                                      "profile/humann2/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        log:
            os.path.join(
                config["output"]["profiling"],
                         "logs/{sample}.humann2.log")
        run:
            if IS_PE:
                shell('''cat {input.reads} > %s''' % input.reads)
                reads = output.reads_merged
            else:
                reads = input.reads[0]

            shell(
                '''
                humann2 \
                --threads {threads} \
                --input %s \
                --output {params.output_dir} \
                --output-basename {params.base_name} \
                --threads {threads} \
                {params.remove} \
                --o-log {log} %s
                ''' % reads)


    rule profiling_humann2_postprocess:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{{sample}}/{{sample}}_{target}.tsv"),
                   target=["genefamilies", "pathabundance", "pathcoverage"])
        output:
            targets = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{{sample}}/{{sample}}_{target}.{norm}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann2"]["normalize_method"]),
            groupprofiles = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{{sample}}/{{sample}}_{group}_groupped.tsv"),
                group=config["params"]["profiling"]["humann2"]["map_database"])

        params:
            normalize_method = config["params"]["profiling"]["humann2"]["normalize_method"],
            regroup_method = config["params"]["profiling"]["humann2"]["regroup_method"],
            map_database =  config["params"]["profiling"]["humann2"]["map_database"]
        run:
            for i in range(0, len(input)):
                shell(
                    '''
                    humann2_renorm_table \
                    --input {input[i]} \
                    --update-snames \
                    --output {output.targets[i]} \
                    --units {params.normalize_method}
                    ''')

            i = -1
            for db in params.map_database:
                i += 1
                shell(
                    '''
                    humann2_regroup_table \
                    --input {input[0]} \
                    --groups %s \
                    --function {params.regroup_method} \
                    --output {output.groupprofiles[i]}
                    '''% db)


    rule profiling_humann2_join:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2/{sample}/{sample}_{target}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2/{sample}/{sample}_{group}_groupped.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   group=config["params"]["profiling"]["humann2"]["map_database"],
                   sample=SAMPLES.index.unique())
        output:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{group}_joined.tsv"),
                group=config["params"]["profiling"]["humann2"]["map_database"])
        params:
            input_dir = os.path.join(config["output"]["profiling"], "profile/humann2"),
            map_database = config["params"]["profiling"]["humann2"]["map_database"]
        run:
            targets = ["genefamilies", "pathabundance", "pathcoverage"]
            for i in range(0, len(targets)):
                shell(
                    '''
                    humann2_join_tables \
                    --input {params.input_dir} \
                    --output {output.targets[i]} \
                    --file_name %s --search-subdirectories
                    ''' % targets[i])

            i = -1
            for db in params.map_database:
                i += 1
                shell(
                    '''
                    humann2_join_tables \
                    --input {params.input_dir} \
                    --output {output.groupprofile[i]} \
                    --file_name %s_groupped --search-subdirectories
                    ''' % db)


    rule profiling_humann2_split_straified:
        input:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{group}_joined.tsv"),
                group=config["params"]["profiling"]["humann2"]["map_database"])
        output:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/{group}_joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   group=config["params"]["profiling"]["humann2"]["map_database"],
                   suffix=["straified", "unstraified"])
        params:
            output_dir = os.path.join(config["output"]["profiling"], "profile"),
            map_database = config["params"]["profiling"]["humann2"]["map_database"]
        run:
            for i in input.targets:
                shell(
                    '''
                    humann2_split_straified_table \
                    -i %s \
                    -o {params.outdir}
                    ''' % i)

            i = -1
            for db in params.map_database:
                i += 1
                shell(
                    '''
                    humann2_split_straified_table \
                    -i {input.groupprofile[i]} \
                    -o {params.output_dir}
                    ''')


    rule profiling_humann2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{group}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/{group}_joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   group=config["params"]["profiling"]["humann2"]["map_database"],
                   suffix=["straified", "unstraified"])

else:
    rule profiling_humann2_all:
        input:


rule profiling_all:
    input:
        rules.profiling_metaphlan2_all.input,
        rules.profiling_jgi_all.input,
        rules.profiling_humann2_all.input,

        rules.rmhost_all.input,
        rules.qcreport_all.input
