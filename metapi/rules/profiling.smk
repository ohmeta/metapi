if config["params"]["profiling"]["metaphlan"]["do_v2"]:
    rule profiling_metaphlan2:
        input:
            assembly_input
        output:
           protected(os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"))
        conda:
            config["envs"]["bioenv2"]
        log:
            os.path.join(
                config["output"]["profiling"],
                         "logs/{sample}.metaphlan2.log")
        params:
            sample_id = "{sample}",
            input_type = config["params"]["profiling"]["metaphlan"]["input_type"],
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v2"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                   else "",
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            no_map = config["params"]["profiling"]["metaphlan"]["no_map"],
            bowtie2_out = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.bowtie2.bz2"),
            biom = config["params"]["profiling"]["metaphlan"]["biom"],
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
        threads:
            config["params"]["profiling"]["threads"]
        run:
           metapi.metaphlan_init(2)
           metapi.merge_metaphlan_tables(input, threads, output=output[0])


    rule profiling_metaphlan2_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2.merged.abundance.profile.tsv"),

            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_metaphlan2_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v3"]:
    rule profiling_metaphlan3:
        input:
            assembly_input
        output:
           protected(os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"))
        log:
            os.path.join(
                config["output"]["profiling"],
                         "logs/{sample}.metaphlan3.log")
        params:
            sample_id = "{sample}",
            input_type = config["params"]["profiling"]["metaphlan"]["input_type"],
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v3"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                   else "",
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            no_unknown_estimation = "--no_unknown_estimation" \
                if config["params"]["profiling"]["metaphlan"]["no_unknown_estimation"] \
                   else "",
            no_map = config["params"]["profiling"]["metaphlan"]["no_map"],
            bowtie2_out = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.bowtie2.bz2"),
            biom = config["params"]["profiling"]["metaphlan"]["biom"],
            biom_out = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.biom"),
            legacy_output = "--legacy-output" \
                if config["params"]["profiling"]["metaphlan"]["legacy_output"] \
                   else "",
            cami_format_output = "--CAMI_format_output" \
                if config["params"]["profiling"]["metaphlan"]["cami_format_output"] \
                   else ""
        threads:
            config["params"]["profiling"]["threads"]
        run:
            shell(
                '''
                metaphlan \
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
                {params.legacy_output} \
                {params.cami_format_output} \
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


    rule profiling_metaphlan3_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"),
                   sample=SAMPLES.index.unique())
        output:
            os.path.join(
                config["output"]["profiling"],
                         "profile/metaphlan3.merged.abundance.profile.tsv")
        threads:
            config["params"]["profiling"]["threads"]
        run:
           metapi.metaphlan_init(3)
           metapi.merge_metaphlan_tables(input, threads, output=output[0])


    rule profiling_metaphlan3_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3.merged.abundance.profile.tsv"),

            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_metaphlan3_all:
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
                          "strain"]),

            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_jgi_all:
        input:


if config["params"]["classify"]["kraken2"]["do"] and \
   config["params"]["profiling"]["bracken"]["do"]:
    rule profiling_bracken:
        input:
            os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.report.gz")
        output:
            profile = protected(os.path.join(
                config["output"]["profiling"],
                "profile/bracken/{sample}/{sample}.bracken.profile")),
            report = protected(os.path.join(
                config["output"]["profiling"],
                "profile/bracken/{sample}/{sample}.bracken.report"))
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/{sample}.bracken.log")
        params:
            database = config["params"]["classify"]["kraken2"]["database"],
            reads_len = config["params"]["profiling"]["bracken"]["reads_len"],
            level = config["params"]["profiling"]["bracken"]["level"]
        threads:
            config["params"]["profiling"]["threads"]
        run:
            import os
            output_dir = os.path.dirname(output.profile)
            report = os.path.join(output_dir, "kreport")

            shell('''pigz -p {threads} -d -k -c {input} > %s''' % report)

            shell(
                '''
                bracken \
                -d {params.database} \
                -i %s \
                -o {output.profile} \
                -w {output.report} \
                -r {params.reads_len} \
                -l {params.level} \
                -t {threads} \
                > {log} 2>&1
                ''' % report)

            shell('''rm -rf %s''' % report)


    rule profiling_bracken_merge:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/bracken/{sample}/{sample}.bracken.profile"),
                 sample=SAMPLES.index.unique())
        output:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/bracken.merged.abundance.profile.{level}.tsv"),
                level=config["params"]["profiling"]["bracken"]["level"])
        log:
            os.path.join(config["output"]["profiling"], "logs/bracken.merged.log")
        run:
            shell(
                '''
                combine_bracken_outputs.py \
                --files {input} \
                --names %s \
                --output {output} \
                > {log} 2>&1
                ''' % ",".join(SAMPLES.index.unique()))


    rule profiling_bracken_all:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/bracken.merged.abundance.profile.{level}.tsv"),
                   level=config["params"]["profiling"]["bracken"]["level"]),
           
            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_bracken_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v2"] and \
   config["params"]["profiling"]["humann2"]["do"]:
    rule profiling_humann2_build_chocophlan_pangenome_db:
        input:
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv")
        output:
            expand(os.path.join(
                config["output"]["profiling"],
                "database/humann2/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
                   suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        conda:
            config["envs"]["bioenv2"]
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/{sample}.humann2.build_pandb.log")
        params:
            id = "{sample}",
            db_dir = os.path.join(config["output"]["profiling"], "database/humann2/{sample}"),
            presense_threshold = config["params"]["profiling"]["humann2"]["presense_threshold"]
        run:
            import logging
            import time
            from humann2 import config
            from humann2.humann2 import timestamp_message
            from humann2.search import prescreen
            from humann2.search import nucleotide

            logger=logging.getLogger(__name__)
            logging.basicConfig(filename=log,
                                format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
                                level=getattr(logging, config.log_level),
                                filemode='w',
                                datefmt='%m/%d/%Y %I:%M:%S %p')
            config.file_basename = params.id
            config.temp_dir = params.db_dir
            config.presense_threshold = params.presense_threshold

            start_time = time.time()
            custom_database = prescreen.create_custom_database(config.nucleotide_database,
                                                               input.profile)
            start_time = timestamp_message("parse taxonomy profile and custom database creation",
                                           start_time)

            nucleotide_index_file = nucleotide.index(custom_database)
            start_time = timestamp_message("database index", start_time)

            shell('''rm -rf %s''' % custom_database)


    rule profiling_humann2:
        input:
            reads = assembly_input,
            index = expand(os.path.join(
                config["output"]["profiling"],
                "database/humann2/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            genefamilies = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_genefamilies.tsv")),
            pathabundance = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_pathabundance.tsv")),
            pathcoverage = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{sample}/{sample}_pathcoverage.tsv"))
        conda:
            config["envs"]["bioenv2"]
        params:
            basename = "{sample}",
            index = os.path.join(config["output"]["profiling"],
                                 "database/humann2/{sample}/{sample}_bowtie2_index"),
            evalue = config["params"]["profiling"]["humann2"]["evalue"],
            presense_threshold = config["params"]["profiling"]["humann2"]["presense_threshold"],
            identity_threshold = config["params"]["profiling"]["humann2"]["identity_threshold"],
            translated_subject_coverage_threshold = \
                config["params"]["profiling"]["humann2"]["translated_subject_coverage_threshold"],
            translated_query_coverage_threshold = \
                config["params"]["profiling"]["humann2"]["translated_query_coverage_threshold"],
            remove = "--remove-temp-output" \
                if config["params"]["profiling"]["humann2"]["remove_temp_output"] \
                   else "",
            memory_use = config["params"]["profiling"]["humann2"]["memory_use"],
            output_dir = os.path.join(config["output"]["profiling"],
                                      "profile/humann2/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        log:
            os.path.join(
                config["output"]["profiling"],
                "logs/{sample}.humann2.log")
        run:
            shell(
                '''
                zcat %s |
                sed 's/ //g' |
                bowtie2 \
                --threads {threads} \
                -x {params.index} \
                -U - | \
                humann2 \
                --threads {threads} \
                --input - \
                --input-format sam \
                --evalue {params.evalue} \
                --presense-threshold {params.presense_threshold} \
                --identity-threshold {params.identity_threshold} \
                --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
                --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
                --memory-use {params.memory_use} \
                --output-basename {params.basename} \
                --output {params.output_dir} \
                {params.remove} \
                --o-log {log} %s
                ''' % " ".join(input.reads))


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
        conda:
            config["envs"]["bioenv2"]
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
        conda:
            config["envs"]["bioenv2"]
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
        conda:
            config["envs"]["bioenv2"]
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
                   suffix=["straified", "unstraified"]),

            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_humann2_all:
        input:


rule profiling_all:
    input:
        rules.profiling_metaphlan2_all.input,
        rules.profiling_metaphlan3_all.input,
        rules.profiling_jgi_all.input,
        rules.profiling_bracken_all.input,
        rules.profiling_humann2_all.input
