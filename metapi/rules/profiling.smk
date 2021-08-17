if config["params"]["profiling"]["bgi_soap"]["do"]:
    rule profiling_bgi_soap:
        input:
            reads = assembly_input_with_short_reads,
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["profiling"]["bgi_soap"]["index_prefix"],
                           suffix=["amb", "ann", "bwt", "fmv", "hot", "lkt", "pac",
                                   "rev.bwt", "rev.fmv", "rev.lkt", "rev.pac", "sa", ".sai"])
        output:
            soap = os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap/{sample}/{sample}.bgi_soap.soap.gz")
        log:
            os.path.join(config["output"]["profiling"], "logs/{sample}.bgi_soap.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/bgi_soap/{sample}.profiling.bgi_soap.benchmark.txt")
        params:
            index = config["params"]["profiling"]["bgi_soap"]["index_prefix"],
            minimal_insert_size = config["params"]["profiling"]["bgi_soap"]["minimal_insert_size"],
            maximal_insert_size = config["params"]["profiling"]["bgi_soap"]["maximal_insert_size"],
            report_repeat_hits = config["params"]["profiling"]["bgi_soap"]["report_repeat_hits"],
            match_model = config["params"]["profiling"]["bgi_soap"]["match_model"],
            align_seed = config["params"]["profiling"]["bgi_soap"]["align_seed"],
            max_mismatch_num = config["params"]["profiling"]["bgi_soap"]["max_mismatch_num"],
            identity = config["params"]["profiling"]["bgi_soap"]["identity"],
            soap = os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap/{sample}/{sample}.bgi_soap.soap")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            if IS_PE:
                shell(
                    '''
                    soap2.22 -a {input.reads[0]} -b {input.reads[1]} -D {params.index} \
                    -m {params.minimal_insert_size} \
                    -x {params.maximal_insert_size} \
                    -r {params.report_repeat_hits} \
                    -l {params.align_seed} \
                    -M {params.match_model} \
                    -v {params.max_mismatch_num} \
                    -c {params.identity} \
                    -S -p {threads} \
                    -o {params.soap}.pe \
                    -2 {params.soap}.se \
                    2> {log}
                    ''')
                shell(
                    '''
                    cat {params.soap}.pe {params.soap}.se > {params.soap}
                    rm -rf {params.soap}.pe {params.soap}.se
                    ''')
            else:
                shell(
                    '''
                    soap2.22 -a {input.reads[0]} -D {params.index} \
                    -m {params.minimal_insert_size} \
                    -x {params.maximal_insert_size} \
                    -r {params.report_repeat_hits} \
                    -l {params.align_seed} \
                    -M {params.match_model} \
                    -v {params.max_mismatch_num} \
                    -c {params.identity} \
                    -S -p {threads} \
                    -o {params.soap} \
                    2> {log}
                    ''')
            shell('''pigz {params.soap}''')


    rule profiling_bgi_soap_merge:
        input:
            soap_list = expand(os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap/{sample}/{sample}.bgi_soap.soap.gz"),
                   sample=SAMPLES.index.unique()),
            taxonomy = config["params"]["profiling"]["index_taxonomy"]
        output:
            abun_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap.merged.abundance.profile.tsv"),
            count_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap.merged.count.profile.tsv")
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/profiling_bgi_soap_merge.log")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            metapi.profiler_init(input.taxonomy)

            count_df, abun_df = metapi.get_all_abun_df(input.soap_list, threads, "bgi_soap")
            count_df.to_csv(output.count_profile, sep='\t', index=False)
            abun_df.to_csv(output.abun_profile, sep='\t', index=False)


    rule profiling_bgi_soap_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap.merged.abundance.profile.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/bgi_soap.merged.count.profile.tsv"),

            rules.qcreport_all.input

else:
    rule profiling_bgi_soap_all:
        input:


if config["params"]["profiling"]["bowtie2"]["do"]:
    rule profiling_bowtie2:
        input:
            reads = assembly_input_with_short_reads,
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["profiling"]["bowtie2"]["index_prefix"],
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]),
            taxonomy = config["params"]["profiling"]["index_taxonomy"]
        output:
            flagstat = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.flagstat"),
            bam = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.sorted.bam"),
            bai = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.sorted.bam.bai"),
            count_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.count.tsv"),
            abun_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.abun.tsv")
        log:
            os.path.join(config["output"]["profiling"], "logs/{sample}.bowtie2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/bowtie2/{sample}.profiling.bowtie2.benchmark.txt")
        params:
            index_prefix = config["params"]["profiling"]["bowtie2"]["index_prefix"]
        threads:
            config["params"]["profiling"]["threads"]
        run:
            shell(
                '''
                bowtie2 \
                -x {params.index_prefix} \
                --end-to-end \
                --very-sensitive \
                --phred33 \
                --threads {threads} \
                --seed 0 \
                --time \
                -k 2 \
                --no-unal \
                --no-discordant \
                %s \
                2> {log} | \
                tee >(samtools flagstat \
                      -@{threads} - \
                      > {output.flagstat}) | \
                samtools sort \
                -m 3G \
                -@{threads} \
                -T {output.bam} \
                -O BAM -o {output.bam} -
                ''' % \
                "-1 {input.reads[0]} -2 {input.reads[1]}" if IS_PE \
                else "-U {input.reads[0]}")

            shell(
                '''
                samtools index -@{threads} {output.bam} {output.bai} 2>> {log}
                ''')

            metapi.profiler_init(input.taxonomy)
            count_df, abun_df = metapi.get_abun_df_bowtie2(output.bam)
            count_df.reset_index().to_csv(output.count_profile, sep='\t', index=False)
            abun_df.reset_index().to_csv(output.abun_profile, sep='\t', index=False)


    rule profiling_bowtie2_merge:
        input:
            count_tsv_list = expand(os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.count.tsv"),
                                    sample=SAMPLES.index.unique()),
            abun_tsv_list = expand(os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2/{sample}/{sample}.bowtie2.abun.tsv"),
                                   sample=SAMPLES.index.unique())
        output:
            abun_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2.merged.abundance.profile.tsv"),
            count_profile = os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2.merged.count.profile.tsv")
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/profiling_bowtie2_merge.log")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            def parse(tsv_file):
                df = pd.read_csv(tsv_file, sep='\t').set_index("lineages_full")
                return df

            count_list = [parse(i) for i in input.count_tsv_list]
            abun_list = [parse(i) for i in input.abun_tsv_list]

            count_df = pd.concat(count_list, axis=1).reset_index()
            abun_df = pd.concat(abun_list, axis=1).reset_index()

            count_df.to_csv(output.count_profile, sep='\t', index=False)
            abun_df.to_csv(output.abun_profile, sep='\t', index=False)


    rule profiling_bowtie2_all:
        input:
            os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2.merged.abundance.profile.tsv"),
            os.path.join(
                config["output"]["profiling"],
                "profile/bowtie2.merged.count.profile.tsv"),

            rules.qcreport_all.input

else:
    rule profiling_bowtie2_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v2"]:
    rule profiling_metaphlan2:
        input:
            reads = assembly_input_with_short_reads
        output:
           profile = protected(os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"))
        conda:
            config["envs"]["bioenv2"]
        log:
            os.path.join(config["output"]["profiling"], "logs/{sample}.metaphlan2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/metaphlan2/{sample}.metaphlan2.benchmark.txt")
        params:
            sample_id = "{sample}",
            wrapper_dir = WRAPPER_DIR,
            input_str = lambda wildcards: ",".join(assembly_input_with_short_reads(wildcards)),
            read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
            bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
            index = config["params"]["profiling"]["metaphlan"]["index_v2"],
            bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
            min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
            taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
            avoid_disqm = "--avoid_disqm" \
                if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                   else "",
            stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v2"],
            stat = config["params"]["profiling"]["metaphlan"]["stat"],
            analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
            bowtie2out = os.path.join(
               config["output"]["profiling"],
               "profile/metaphlan2/{sample}/{sample}.metaphlan2.bowtie2.bz2"),
            save_bowtie2out = config["params"]["profiling"]["metaphlan"]["save_bowtie2out"]
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            python {params.wrapper_dir}/metaphlan2_wrapper.py \
            --input {input.reads} \
            --analysis_type {params.analysis_type} \
            --input_type multifastq \
            --bowtie2db {params.bowtie2db} \
            --index {params.index} \
            --bt2_ps {params.bowtie2_presets} \
            --bowtie2out {params.bowtie2out} \
            --tax_lev {params.taxonomic_level} \
            --min_cu_len {params.min_cu_len} \
            --stat_q {params.stat_q} \
            --output_file {output.profile} \
            --sample_id {params.sample_id} \
            --nproc {threads} \
            --read_min_len {params.read_min_len} \
            >{log} 2>&1
            '''


    rule profiling_metaphlan2_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan2/{sample}/{sample}.metaphlan2.abundance.profile.tsv"),
                   sample=SAMPLES.index.unique())
        output:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan2.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"])
        threads:
            config["params"]["profiling"]["threads"]
        run:
           metapi.metaphlan_init(2)
           df_list = metapi.merge_metaphlan_tables(input, threads)
           for i in range(0, len(df_list)):
               df_list[i].to_csv(output[i], sep='\t', index=False)


    rule profiling_metaphlan2_all:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan2.merged.abundance.profile.{level}.tsv"),
                level=["all", "superkingdom", "phylum", "class",
                       "order", "family", "genus", "species", "strain"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_metaphlan2_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v3"]:
    if not config["params"]["profiling"]["metaphlan"]["do_v3_one_way"]:
        rule profiling_metaphlan3:
            input:
                reads = assembly_input_with_short_reads,
                index_mp3 = expand(
                    os.path.join(
                        config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                        "{mpa_name}.{suffix}"),
                    mpa_name = config["params"]["profiling"]["metaphlan"]["index_v3"],
                    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"])
            output:
               protected(os.path.join(
                   config["output"]["profiling"],
                   "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"))
            log:
                os.path.join(config["output"]["profiling"], "logs/{sample}.metaphlan3.log")
            benchmark:
                os.path.join(config["output"]["profiling"],
                            "benchmark/metaphlan3/{sample}.metaphlan3.benchmark.txt")
            conda:
                config["envs"]["biobakery"]
            params:
                sample_id = "{sample}",
                read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
                bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                index = config["params"]["profiling"]["metaphlan"]["index_v3"],
                bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
                min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
                taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
                avoid_disqm = "--avoid_disqm" \
                    if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                    else "",
                stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v3"],
                stat = config["params"]["profiling"]["metaphlan"]["stat"],
                analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
                unknown_estimation = "--unknown_estimation" \
                    if config["params"]["profiling"]["metaphlan"]["unknown_estimation"] \
                    else "",
                add_viruses = "--add_viruses" \
                    if config["params"]["profiling"]["metaphlan"]["add_viruses"] \
                    else "",
                map_out = "--no_map" if config["params"]["profiling"]["metaphlan"]["no_map"] \
                    else "--bowtie2out %s" % os.path.join(
                        config["output"]["profiling"],
                        "profile/metaphlan3/{sample}/{sample}.metaphlan3.bowtie2.bz2"),
                biom_out = "--biom %s" % os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.biom") \
                    if config["params"]["profiling"]["metaphlan"]["biom"] \
                    else "",
                legacy_output = "--legacy-output" \
                    if config["params"]["profiling"]["metaphlan"]["legacy_output"] \
                    else "",
                cami_format_output = "--CAMI_format_output" \
                    if config["params"]["profiling"]["metaphlan"]["cami_format_output"] \
                    else ""
            threads:
                config["params"]["profiling"]["threads"]
            shell:
                '''
                reads=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.reads})

                metaphlan \
                $reads \
                --input_type fastq \
                --read_min_len {params.read_min_len} \
                --nproc {threads} \
                {params.map_out} \
                --bowtie2db {params.bowtie2db} \
                --index {params.index} \
                --bt2_ps {params.bowtie2_presets} \
                --min_cu_len {params.min_cu_len} \
                --tax_lev {params.taxonomic_level} \
                {params.avoid_disqm} \
                --stat_q {params.stat_q} \
                --stat {params.stat} \
                -t {params.analysis_type} \
                {params.unknown_estimation} \
                {params.add_viruses} \
                --output_file {output} \
                --sample_id {params.sample_id} \
                {params.legacy_output} \
                {params.cami_format_output} \
                {params.biom_out} \
                2> {log}
                '''
    else:
        rule profiling_metaphlan3:
            input:
                reads = lambda wildcards: get_reads(wildcards, "raw"),
                index_mp3 = expand(
                    os.path.join(
                        config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                        "{mpa_name}.{suffix}"),
                    mpa_name = config["params"]["profiling"]["metaphlan"]["index_v3"],
                    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "pkl"]),
                index_host = expand("{prefix}.{suffix}",
                    prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                    suffix=BWA_INDEX_SUFFIX)
            output:
                profile = protected(os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")),
                html = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.fastp.html"),
                json = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.fastp.json"),
                flagstat_rmhost = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.rmhost.flagstat")
            log:
                os.path.join(config["output"]["profiling"],
                             "logs/metaphlan3/{sample}.metaphlan3.log")
            benchmark:
                os.path.join(config["output"]["profiling"],
                            "benchmark/metaphlan3/{sample}.metaphlan3.benchmark.txt")
            #conda:
            #    config["envs"]["biobakery"]
            params:
                bwa = "bwa-mem2" if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" else "bwa",
                minimum_seed_length = config["params"]["rmhost"]["bwa"]["minimum_seed_length"],
                host_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
                trimming_log = os.path.join(config["output"]["profiling"],
                                            "logs/fastp/{sample}.fastp.log"),
                rmhost_log = os.path.join(config["output"]["profiling"],
                                          "logs/bwa/{sample}.bwa.log"),
                sample_id = "{sample}",
                read_min_len = config["params"]["profiling"]["metaphlan"]["read_min_len"],
                bowtie2db = config["params"]["profiling"]["metaphlan"]["bowtie2db"],
                index = config["params"]["profiling"]["metaphlan"]["index_v3"],
                bowtie2_presets = config["params"]["profiling"]["metaphlan"]["bowtie2_presets"],
                min_cu_len = config["params"]["profiling"]["metaphlan"]["min_cu_len"],
                taxonomic_level = config["params"]["profiling"]["metaphlan"]["taxonomic_level"],
                avoid_disqm = "--avoid_disqm" \
                    if config["params"]["profiling"]["metaphlan"]["avoid_disqm"] \
                    else "",
                stat_q = config["params"]["profiling"]["metaphlan"]["stat_q_v3"],
                stat = config["params"]["profiling"]["metaphlan"]["stat"],
                analysis_type = config["params"]["profiling"]["metaphlan"]["analysis_type"],
                unknown_estimation = "--unknown_estimation" \
                    if config["params"]["profiling"]["metaphlan"]["unknown_estimation"] \
                    else "",
                map_out = "--no_map" if config["params"]["profiling"]["metaphlan"]["no_map"] \
                    else "--bowtie2out %s" % os.path.join(
                        config["output"]["profiling"],
                        "profile/metaphlan3/{sample}/{sample}.metaphlan3.bowtie2.bz2"),
                biom_out = "--biom %s" % os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.biom") \
                    if config["params"]["profiling"]["metaphlan"]["biom"] \
                    else "",
                legacy_output = "--legacy-output" \
                    if config["params"]["profiling"]["metaphlan"]["legacy_output"] \
                    else "",
                cami_format_output = "--CAMI_format_output" \
                    if config["params"]["profiling"]["metaphlan"]["cami_format_output"] \
                    else ""
            threads:
                config["params"]["profiling"]["threads"]
            run:
                shell(f'''mkdir -p {os.path.dirname(params.trimming_log)}''')
                shell(f'''mkdir -p {os.path.dirname(params.rmhost_log)}''')
                if TRIMMING_DO and RMHOST_DO:
                    shell('''date > {log}''')
                    shell(
                        '''
                        fastp %s %s --thread {threads} \
                        --stdout --json {output.json} --html {output.html} 2> {params.trimming_log} | \
                        {params.bwa} mem -p -t {threads} \
                        -k {params.minimum_seed_length} {params.host_prefix} - 2> {params.rmhost_log} | \
                            tee >(samtools flagstat \
                                  -@{threads} - \
                                  > {output.flagstat_rmhost}) | \
                        samtools fastq -@{threads} -N -f 12 -F 256 - | \
                        metaphlan \
                        --input_type fastq \
                        --read_min_len {params.read_min_len} \
                        --nproc {threads} \
                        {params.map_out} \
                        --bowtie2db {params.bowtie2db} \
                        --index {params.index} \
                        --bt2_ps {params.bowtie2_presets} \
                        --min_cu_len {params.min_cu_len} \
                        --tax_lev {params.taxonomic_level} \
                        {params.avoid_disqm} \
                        --stat_q {params.stat_q} \
                        --stat {params.stat} \
                        -t {params.analysis_type} \
                        {params.unknown_estimation} \
                        --output_file {output.profile} \
                        --sample_id {params.sample_id} \
                        {params.legacy_output} \
                        {params.cami_format_output} \
                        {params.biom_out} \
                        2>> {log}
                        ''' % (
                            "--in1 {input.reads[0]} --in2 {input.reads[1]}" if IS_PE else "--in1 {input.reads[0]}",
                            f"{ADAPTER_OPERATION}"
                        ))
                    shell('''date >> {log}''')


        rule profiling_metaphlan3_oneway_qc:
            input:
                jsons = expand(os.path.join(
                    config["output"]["profiling"],
                    "stats_preprocess/{sample}/{sample}.fastp.json"),
                    sample=SAMPLES.index.unique()),
                flagstats = expand(os.path.join(
                    config["output"]["profiling"],
                    "stats_preprocess/{sample}/{sample}.rmhost.flagstat"),
                    sample=SAMPLES.index.unique())
            output:
                qc = os.path.join(config["output"]["profiling"], "qc_stats.tsv"),
            run:
                import concurrent.futures

                fastp_dict_list = []
                json_tuples = [(j, IS_PE) for j in input.jsons]
                with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
                    for fastp_dict in executor.map(metapi.parse_fastp_json, json_tuples):
                        fastp_dict_list.append(fastp_dict)
                fastp_df = pd.DataFrame(fastp_dict_list)
                flagstat_df = metapi.flagstats_summary(input.flagstats, 2)
                qc_df = fastp_df.merge(faststat_df)
                qc_df["after_rmhosted_total_reads"] = qc_df["total_num"] - qc_df["mapped_num"] - qc_df["singletons_num"]
                qc_df["host_rate"] = (qc_df["mapped_num"] + qc_df["singletons_num"]) / qc_df["total_num"]

                qc_df.loc[:, list(fastp_df.columns) + ["after_rmhosted_total_reads", "host_rate"]]\
                    .to_csv(output.qc, sep="\t", index=False)


    rule profiling_metaphlan3_merge:
        input:
           abuns = expand(os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv"),
                          sample=SAMPLES.index.unique())
        output:
            profiles = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/metaphlan3.merged.abundance.profile.{level}.tsv"),
                level=["all", "strain", "species", "genus", "family",
                       "order", "class", "phylum", "superkingdom"])
        threads:
            config["params"]["profiling"]["threads"]
        run:
            metapi.metaphlan_init(3)
            profile_list = metapi.merge_metaphlan_tables(input.abuns, threads)
            for i in range(0, len(profile_list)):
                profile_list[i].to_csv(output.profiles[i], sep='\t', index=False)


    if config["params"]["profiling"]["metaphlan"]["do_v3_one_way"]:
        rule profiling_metaphlan3_all:
            input:
                os.path.join(config["output"]["profiling"], "qc_stats.tsv"),
                expand(
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/metaphlan3.merged.abundance.profile.{level}.tsv"),
                    level=["all", "superkingdom", "phylum", "class",
                           "order", "family", "genus", "species", "strain"])

    else:
        rule profiling_metaphlan3_all:
            input:
                expand(
                    os.path.join(
                        config["output"]["profiling"],
                        "profile/metaphlan3.merged.abundance.profile.{level}.tsv"),
                    level=["all", "superkingdom", "phylum", "class",
                           "order", "family", "genus", "species", "strain"]),
                
                #rules.rmhost_all.input,
                rules.qcreport_all.input

else:
    rule profiling_metaphlan3_all:
        input:
          

if config["params"]["profiling"]["jgi"]["do"]:
    if not config["params"]["profiling"]["jgi"]["oneway"]:
        rule profiling_jgi:
            input:
                reads = assembly_input_with_short_reads,
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
                os.path.join(config["output"]["profiling"], "logs/{sample}.jgi.log")
            benchmark:
                os.path.join(config["output"]["profiling"],
                             "benchmark/jgi/{sample}.jgi.benchmark.txt")
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
                coverage = protected(os.path.join(
                    config["output"]["profiling"],
                    "profile/jgi/{sample}/{sample}.jgi.coverage.gz")),
                html = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.fastp.html"),
                json = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.fastp.json"),
                flagstat_rmhost = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.rmhost.flagstat"),
                flagstat_profiling = os.path.join(config["output"]["profiling"],
                                    "stats_preprocess/{sample}/{sample}.profiling.flagstat")
            log:
                os.path.join(config["output"]["profiling"], "logs/{sample}.trimming.rmhost.jgi.log")
            benchmark:
                os.path.join(config["output"]["profiling"],
                             "benchmark/jgi_2/{sample}.jgi.benchmark.txt")
            threads:
                config["params"]["profiling"]["threads"]
            params:
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
                        fastp %s %s --thread {threads} \
                        --stdout --json {output.json} --html {output.html} 2>> {log} | \
                        bowtie2 --threads {threads} -x {params.host_prefix} %s - 2>> {log} | \
                            tee >(samtools flagstat \
                                  -@{threads} - \
                                  > {output.flagstat_rmhost}) | \
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
                            tee >(samtools flagstat \
                                  -@{threads} - \
                                  > {output.flagstat_profiling}) | \
                        sambamba view -q --nthreads {threads} \
                        --compression-level {params.compression_level} \
                        --format bam \
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
                            f"{ADAPTER_OPERATION}",
                            "--interleaved" if IS_PE else "",
                            "--interleaved" if IS_PE else ""
                        ))
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

    rule profiling_jgi_all_:
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


    if not config["params"]["profiling"]["jgi"]["oneway"]:
        rule profiling_jgi_all:
            input:
                rules.profiling_jgi_all_.input,

                #rules.rmhost_all.input,
                rules.qcreport_all.input

    else:
        rule profiling_jgi_all:
            input:
                rules.profiling_jgi_all_.input

else:
    rule profiling_jgi_all:
        input:


if config["params"]["classify"]["kraken2"]["do"] and \
   config["params"]["profiling"]["bracken"]["do"]:
    rule profiling_bracken:
        input:
            os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.report")
        output:
            profile = protected(os.path.join(
                config["output"]["profiling"],
                "profile/bracken/{sample}/{sample}.bracken.{level}.profile")),
            report = protected(os.path.join(
                config["output"]["profiling"],
                "profile/bracken/{sample}/{sample}.bracken.{level}.report"))
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/{sample}.bracken.{level}.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/bracken/{sample}.bracken.{level}.benchmark.txt")
        params:
            database = config["params"]["classify"]["kraken2"]["database"],
            reads_len = config["params"]["profiling"]["bracken"]["reads_len"],
            level = "{level}"
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            bracken \
            -d {params.database} \
            -i {input} \
            -o {output.profile} \
            -w {output.report} \
            -r {params.reads_len} \
            -l {params.level} \
            -t {threads} \
            > {log} 2>&1
            '''


    rule profiling_bracken_merge:
        input:
            expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/bracken/{sample}/{sample}.bracken.{{level}}.profile"),
                 sample=SAMPLES.index.unique())
        output:
            os.path.join(
                config["output"]["profiling"],
                "profile/bracken.merged.abundance.profile.{level}.tsv")
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/bracken.merged.{level}.log")
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
           
            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_bracken_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v2"] and \
   config["params"]["profiling"]["humann"]["do_v2"]:
    if config["params"]["profiling"]["humann"]["update_config"]:
        rule profiling_humann2_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann2.config.done"))
            log:
                os.path.join(config["output"]["profiling"], "logs/humann2.config.log")
            conda:
                config["envs"]["bioenv2"]
            params:
                database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping"],
                database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide"],
                database_protein = config["params"]["profiling"]["humann"]["database_protein"],
                threads = config["params"]["profiling"]["threads"]
            shell:
                '''
                humann2_config > {log}
                humann2_config --update database_folders utility_mapping {params.database_utility_mapping}
                humann2_config --update database_folders nucleotide {params.database_nucleotide}
                humann2_config --update database_folders protein {params.database_protein}
                humann2_config --update run_modes threads {params.threads}
                echo "####" >> {log}
                humann2_config >> {log}
                '''
    else:
        rule profiling_humann2_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann2.config.done"))
            shell:
                '''
                echo "hello"
                '''


    rule profiling_humann2_build_chocophlan_pangenome_db:
        input:
            tag = os.path.join(config["output"]["profiling"], ".humann2.config.done"),
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
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/humann2/{sample}.bowtie2_index.benchmark.txt")
        params:
            basename = "{sample}",
            wrapper_dir = WRAPPER_DIR,
            db_dir = os.path.join(config["output"]["profiling"], "database/humann2/{sample}"),
            prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"]
        shell:
            '''
            python {params.wrapper_dir}/humann2_db_wrapper.py \
            --log {log} \
            --basename {params.basename} \
            --db_dir {params.db_dir} \
            --prescreen_threshold {params.prescreen_threshold} \
            --taxonomic_profile {input.profile}
            '''


    rule profiling_humann2:
        input:
            tag = os.path.join(config["output"]["profiling"], ".humann2.config.done"),
            reads = assembly_input_with_short_reads,
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
        log:
            os.path.join(config["output"]["profiling"], "logs/{sample}.humann2.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/humann2/{sample}.humann2.benchmark.txt")
        conda:
            config["envs"]["bioenv2"]
        params:
            basename = "{sample}",
            index = os.path.join(config["output"]["profiling"],
                                 "database/humann2/{sample}/{sample}_bowtie2_index"),
            evalue = config["params"]["profiling"]["humann"]["evalue"],
            prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
            identity_threshold = config["params"]["profiling"]["humann"]["identity_threshold"],
            translated_subject_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
            translated_query_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
            xipe = config["params"]["profiling"]["humann"]["xipe"],
            minpath = config["params"]["profiling"]["humann"]["minpath"],
            pick_frames = config["params"]["profiling"]["humann"]["pick_frames"],
            gap_fill = config["params"]["profiling"]["humann"]["gap_fill"],
            remove_temp_output = "--remove-temp-output" \
                if config["params"]["profiling"]["humann"]["remove_temp_output"] \
                   else "",
            memory_use = config["params"]["profiling"]["humann"]["memory_use"],
            output_dir = os.path.join(config["output"]["profiling"],
                                      "profile/humann2/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            zcat {input.reads} | \
            bowtie2 \
            --threads {threads} \
            -x {params.index} \
            -U - 2>> {log} | \
            humann2 \
            --threads {threads} \
            --input - \
            --input-format sam \
            --evalue {params.evalue} \
            --prescreen-threshold {params.prescreen_threshold} \
            --identity-threshold {params.identity_threshold} \
            --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
            --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
            --xipe {params.xipe} \
            --minpath {params.minpath} \
            --pick-frames {params.pick_frames} \
            --gap-fill {params.gap_fill} \
            --memory-use {params.memory_use} \
            --output-basename {params.basename} \
            --output {params.output_dir} \
            {params.remove_temp_output} \
            --o-log {log}
            '''


    rule profiling_humann2_postprocess:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{{sample}}/{{sample}}_{target}.tsv"),
                   target=["genefamilies", "pathabundance", "pathcoverage"])
        output:
            targets = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"]),
            groupprofiles = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann2/{{sample}}/{{sample}}_{group}_groupped.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/{sample}.humann2_postprocess.log")
        conda:
            config["envs"]["bioenv2"]
        params:
            wrapper_dir =WRAPPER_DIR,
            normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
            regroup_method = config["params"]["profiling"]["humann"]["regroup_method"],
            map_database =  config["params"]["profiling"]["humann"]["map_database"]
        shell:
            '''
            humann2_renorm_table \
            --input {input[0]} \
            --update-snames \
            --output {output.targets[0]} \
            --units {params.normalize_method} \
            > {log} 2>&1

            humann2_renorm_table \
            --input {input[1]} \
            --update-snames \
            --output {output.targets[1]} \
            --units {params.normalize_method} \
            >> {log} 2>&1

            humann2_renorm_table \
            --input {input[2]} \
            --update-snames \
            --output {output.targets[2]} \
            --units {params.normalize_method} \
            >> {log} 2>&1

            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            regroup_table \
            --input {input[0]} \
            --groups {params.map_database} \
            --function {params.regroup_method} \
            --output {output.groupprofiles} \
            >> {log} 2>&1
            '''


    rule profiling_humann2_join:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2/{sample}/{sample}_{target}.tsv"),
               os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2/{sample}/{sample}_{target}_{norm}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2/{sample}/{sample}_{group}_groupped.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   sample=SAMPLES.index.unique())
        output:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            targets_norm = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_{norm}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2.{group}.joined.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann2_join.log")
        conda:
            config["envs"]["bioenv2"]
        params:
            wrapper_dir =WRAPPER_DIR,
            input_dir = os.path.join(config["output"]["profiling"], "profile/humann2"),
            normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
            map_database = config["params"]["profiling"]["humann"]["map_database"]
        shell:
            '''
            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            join_tables \
            --input {params.input_dir} \
            --output {output.targets} \
            --file_name genefamilies.tsv pathabundance.tsv pathcoverage.tsv \
            > {log} 2>&1

            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            join_tables \
            --input {params.input_dir} \
            --output {output.targets_norm} \
            --file_name \
            genefamilies_{params.normalize_method}.tsv \
            pathabundance_{params.normalize_method}.tsv \
            pathcoverage_{params.normalize_method}.tsv \
            > {log} 2>&1

            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            join_tables \
            --input {params.input_dir} \
            --output {output.groupprofile} \
            --file_name {params.map_database} \
            >> {log} 2>&1
            '''


    rule profiling_humann2_split_stratified:
        input:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            targets_norm = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_{norm}_joined.tsv"),
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2.{group}.joined.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        output:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann2_split_stratified.log")
        conda:
            config["envs"]["bioenv2"]
        params:
            wrapper_dir = WRAPPER_DIR,
            output_dir = os.path.join(config["output"]["profiling"], "profile"),
            map_database = config["params"]["profiling"]["humann"]["map_database"]
        shell:
            '''
            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            split_stratified_table \
            --input {input.targets} \
            --output {params.output_dir} \
            > {log} 2>&1

            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            split_stratified_table \
            --input {input.targets_norm} \
            --output {params.output_dir} \
            > {log} 2>&1

            python {params.wrapper_dir}/humann2_postprocess_wrapper.py \
            split_stratified_table \
            --input {input.groupprofile} \
            --output {params.output_dir} \
            >> {log} 2>&1
            '''


    rule profiling_humann2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_joined_{suffix}.tsv"),
               os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_{norm}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2.{group}.joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann2.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_humann2_all:
        input:


if config["params"]["profiling"]["metaphlan"]["do_v3"] and \
   config["params"]["profiling"]["humann"]["do_v3"]:
    if config["params"]["profiling"]["humann"]["update_config"]:
        rule profiling_humann3_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann3.config.done"))
            log:
                os.path.join(config["output"]["profiling"], "logs/humann3.config.log")
            conda:
                config["envs"]["biobakery"]
            params:
                database_utility_mapping = config["params"]["profiling"]["humann"]["database_utility_mapping_v3"],
                database_nucleotide = config["params"]["profiling"]["humann"]["database_nucleotide_v3"],
                database_protein = config["params"]["profiling"]["humann"]["database_protein_v3"],
                threads = config["params"]["profiling"]["threads"]
            shell:
                '''
                humann_config > {log}
                humann_config --update database_folders utility_mapping {params.database_utility_mapping}
                humann_config --update database_folders nucleotide {params.database_nucleotide}
                humann_config --update database_folders protein {params.database_protein}
                humann_config --update run_modes threads {params.threads}
                echo "####" >> {log}
                humann_config >> {log}
                '''
    else:
        rule profiling_humann3_config:
            output:
                touch(os.path.join(config["output"]["profiling"], ".humann3.config.done"))
            shell:
                '''
                echo "hello"
                '''


    #rule profiling_humann3_build_chocophlan_pangenome_db:
    #    input:
    #        profile = os.path.join(
    #            config["output"]["profiling"],
    #            "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")
    #    output:
    #        expand(os.path.join(
    #            config["output"]["profiling"],
    #            "database/humann3/{{sample}}/{{sample}}_bowtie2_index.{suffix}"),
    #               suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    #    log:
    #        os.path.join(config["output"]["profiling"],
    #                     "logs/{sample}.humann3.build_pandb.log")
    #    params:
    #        basename = "{sample}",
    #        wrapper_dir = WRAPPER_DIR,
    #        db_dir = os.path.join(config["output"]["profiling"], "database/humann3/{sample}"),
    #        prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"]
    #    shell:
    #        '''
    #        python {params.wrapper_dir}/humann3_db_wrapper.py \
    #        --log {log} \
    #        --basename {params.basename} \
    #        --db_dir {params.db_dir} \
    #        --prescreen_threshold {params.prescreen_threshold} \
    #        --taxonomic_profile {input.profile}
    #        '''


    rule profiling_humann3:
        input:
            tag = os.path.join(config["output"]["profiling"], ".humann3.config.done"),
            reads = assembly_input_with_short_reads,
            profile = os.path.join(
                config["output"]["profiling"],
                "profile/metaphlan3/{sample}/{sample}.metaphlan3.abundance.profile.tsv")
        output:
            genefamilies = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{sample}/{sample}_genefamilies.tsv")),
            pathabundance = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{sample}/{sample}_pathabundance.tsv")),
            pathcoverage = protected(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{sample}/{sample}_pathcoverage.tsv"))
        log:
            os.path.join(config["output"]["profiling"], "logs/{sample}.humann3.log")
        benchmark:
            os.path.join(config["output"]["profiling"],
                         "benchmark/humann3/{sample}.humann3.benchmark.txt")
        conda:
            config["envs"]["biobakery"]
        params:
            basename = "{sample}",
            wrapper_dir = WRAPPER_DIR,
            index = os.path.join(config["output"]["profiling"],
                                 "database/humann3/{sample}/{sample}_bowtie2_index"),
            evalue = config["params"]["profiling"]["humann"]["evalue"],
            prescreen_threshold = config["params"]["profiling"]["humann"]["prescreen_threshold"],
            nucleotide_identity_threshold = \
                config["params"]["profiling"]["humann"]["nucleotide_identity_threshold"],
            translated_identity_threshold = \
                config["params"]["profiling"]["humann"]["translated_identity_threshold"],
            translated_subject_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_subject_coverage_threshold"],
            translated_query_coverage_threshold = \
                config["params"]["profiling"]["humann"]["translated_query_coverage_threshold"],
            nucleotide_subject_coverage_threshold = \
                config["params"]["profiling"]["humann"]["nucleotide_subject_coverage_threshold"],
            nucleotide_query_coverage_threshold = \
                config["params"]["profiling"]["humann"]["nucleotide_query_coverage_threshold"],
            minpath = config["params"]["profiling"]["humann"]["minpath"],
            gap_fill = config["params"]["profiling"]["humann"]["gap_fill"],
            xipe = config["params"]["profiling"]["humann"]["xipe"],
            pathways = config["params"]["profiling"]["humann"]["pathways"],
            remove_temp_output = "--remove-temp-output" \
                if config["params"]["profiling"]["humann"]["remove_temp_output"] \
                   else "",
            memory_use = config["params"]["profiling"]["humann"]["memory_use"],
            output_dir = os.path.join(config["output"]["profiling"],
                                      "profile/humann3/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            mkdir -p {params.output_dir}

            python {params.wrapper_dir}/misc.py \
            --basename {params.basename} \
            --input-file {input.reads} \
            --output-dir {params.output_dir}

            rm -rf {params.output_dir}/{params.basename}_humann_temp_*

            humann \
            --threads {threads} \
            --input {params.output_dir}/{params.basename}.fq.gz \
            --input-format fastq.gz \
            --taxonomic-profile {input.profile} \
            --evalue {params.evalue} \
            --prescreen-threshold {params.prescreen_threshold} \
            --nucleotide-identity-threshold {params.nucleotide_identity_threshold} \
            --translated-identity-threshold {params.translated_identity_threshold} \
            --translated-subject-coverage-threshold {params.translated_subject_coverage_threshold} \
            --nucleotide-subject-coverage-threshold {params.nucleotide_subject_coverage_threshold} \
            --translated-query-coverage-threshold {params.translated_query_coverage_threshold} \
            --nucleotide-query-coverage-threshold {params.nucleotide_query_coverage_threshold} \
            --gap-fill {params.gap-fill} \
            --minpath {params.minpath} \
            --xipe {params.xipe} \
            --memory-use {params.memory_use} \
            --output-basename {params.basename} \
            --output {params.output_dir} \
            {params.remove_temp_output} \
            --o-log {log}

            rm -rf {params.output_dir}/{params.basename}.fq.gz
            '''


    rule profiling_humann3_postprocess:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{{sample}}/{{sample}}_{target}.tsv"),
                   target=["genefamilies", "pathabundance", "pathcoverage"])
        output:
            targets = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{{sample}}/{{sample}}_{target}_{norm}.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"]),
            groupprofiles = expand(os.path.join(
                config["output"]["profiling"],
                "profile/humann3/{{sample}}/{{sample}}_{group}_groupped.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/{sample}.humann3_postprocess.log")
        conda:
            config["envs"]["biobakery"]
        params:
            wrapper_dir =WRAPPER_DIR,
            normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
            regroup_method = config["params"]["profiling"]["humann"]["regroup_method"],
            map_database =  config["params"]["profiling"]["humann"]["map_database"]
        shell:
            '''
            humann_renorm_table \
            --input {input[0]} \
            --update-snames \
            --output {output.targets[0]} \
            --units {params.normalize_method} \
            > {log} 2>&1

            humann_renorm_table \
            --input {input[1]} \
            --update-snames \
            --output {output.targets[1]} \
            --units {params.normalize_method} \
            >> {log} 2>&1

            humann_renorm_table \
            --input {input[2]} \
            --update-snames \
            --output {output.targets[2]} \
            --units {params.normalize_method} \
            >> {log} 2>&1

            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            regroup_table \
            --input {input[0]} \
            --groups {params.map_database} \
            --function {params.regroup_method} \
            --output {output.groupprofiles} \
            >> {log} 2>&1
            '''


    rule profiling_humann3_join:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3/{sample}/{sample}_{target}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3/{sample}/{sample}_{target}_{norm}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3/{sample}/{sample}_{group}_groupped.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm=config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   sample=SAMPLES.index.unique())
        output:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            targets_norm = expand(
               os.path.join(
                   config["output"]["profiling"],
                   "profile/humann3_{target}_{norm}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"],
                norm=config["params"]["profiling"]["humann"]["normalize_method"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3.{group}.joined.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann3_join.log")
        conda:
            config["envs"]["biobakery"]
        params:
            wrapper_dir =WRAPPER_DIR,
            input_dir = os.path.join(config["output"]["profiling"], "profile/humann3"),
            normalize_method = config["params"]["profiling"]["humann"]["normalize_method"],
            map_database = config["params"]["profiling"]["humann"]["map_database"]
        shell:
            '''
            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            join_tables \
            --input {params.input_dir} \
            --output {output.targets} \
            --file_name genefamilies.tsv pathabundance.tsv pathcoverage.tsv \
            > {log} 2>&1

            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            join_tables \
            --input {params.input_dir} \
            --output {output.targets_norm} \
            --file_name \
            genefamilies_{params.normalize_method}.tsv \
            pathabundance_{params.normalize_method}.tsv \
            pathcoverage_{params.normalize_method}.tsv \
            > {log} 2>&1

            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            join_tables \
            --input {params.input_dir} \
            --output {output.groupprofile} \
            --file_name {params.map_database} \
            >> {log} 2>&1
            '''


    rule profiling_humann3_split_stratified:
        input:
            targets = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_joined.tsv"),
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            targets_norm = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_{norm}_joined.tsv"),
                norm = config["params"]["profiling"]["humann"]["normalize_method"],
                target=["genefamilies", "pathabundance", "pathcoverage"]),
            groupprofile = expand(
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3.{group}.joined.tsv"),
                group=config["params"]["profiling"]["humann"]["map_database"])
        output:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"])
        log:
            os.path.join(config["output"]["profiling"],
                         "logs/humann3_split_stratified.log")
        conda:
            config["envs"]["biobakery"]
        params:
            wrapper_dir = WRAPPER_DIR,
            output_dir = os.path.join(config["output"]["profiling"], "profile"),
            map_database = config["params"]["profiling"]["humann"]["map_database"]
        shell:
            '''
            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            split_stratified_table \
            --input {input.targets} \
            --output {params.output_dir} \
            > {log} 2>&1

            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            split_stratified_table \
            --input {input.targets_norm} \
            --output {params.output_dir} \
            > {log} 2>&1

            python {params.wrapper_dir}/humann3_postprocess_wrapper.py \
            split_stratified_table \
            --input {input.groupprofile} \
            --output {params.output_dir} \
            >> {log} 2>&1
            '''


    rule profiling_humann3_all:
        input:
            expand([
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_{norm}_joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3_{target}_{norm}_joined_{suffix}.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3.{group}.joined.tsv"),
                os.path.join(
                    config["output"]["profiling"],
                    "profile/humann3.{group}.joined_{suffix}.tsv")],
                   target=["genefamilies", "pathabundance", "pathcoverage"],
                   norm = config["params"]["profiling"]["humann"]["normalize_method"],
                   group=config["params"]["profiling"]["humann"]["map_database"],
                   suffix=["stratified", "unstratified"]),

            #rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule profiling_humann3_all:
        input:


rule profiling_all:
    input:
        rules.profiling_bgi_soap_all.input,
        rules.profiling_bowtie2_all.input,
        rules.profiling_metaphlan2_all.input,
        rules.profiling_metaphlan3_all.input,
        rules.profiling_jgi_all.input,
        rules.profiling_bracken_all.input,
        rules.profiling_humann2_all.input,
        rules.profiling_humann3_all.input
