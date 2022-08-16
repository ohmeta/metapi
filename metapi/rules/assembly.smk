def get_reads_for_assembly(wildcards, step, have_single=False, have_long=False):
    samples_id_list = metapi.get_samples_id_by_assembly_group(SAMPLES, wildcards.assembly_group)
    short_reads = get_short_reads_list(step, samples_id_list)

    if have_long:
        long_reads = expand(os.path.join(
            config["output"]["raw"],
            "long_reads/{sample}/{sample}.{step}{read}.fq"),
            step="raw",
            read=long_reads_suffix(),
            sample=samples_id_list)
        return [short_reads, long_reads]
    else:
        return short_reads


def assembly_input_with_short_reads(wildcards):
    if RMHOST_DO:
        return get_reads_for_assembly(wildcards, "rmhost", False, False)
    elif TRIMMING_DO:
        return get_reads_for_assembly(wildcards, "trimming", False, False)
    else:
        return get_reads_for_assembly(wildcards, "raw", False, False)


def assembly_input_with_short_and_long_reads(wildcards):
    if RMHOST_DO:
        return get_reads_for_assembly(wildcards, "rmhost", False, True)
    elif TRIMMING_DO:
        return get_reads_for_assembly(wildcards, "trimming", False, True)
    else:
        return get_reads_for_assembly(wildcards, "raw", False, True)


if "megahit" in ASSEMBLERS:
    rule assembly_megahit:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.gfa.gz"))
        conda:
            config["envs"]["megahit"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/megahit/{assembly_group}.megahit.benchmark.txt")
        priority:
            20
        params:
            output_prefix = "{assembly_group}",
            min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
            k_list = ",".join(config["params"]["assembly"]["megahit"]["k_list"]),
            presets = config["params"]["assembly"]["megahit"]["presets"],
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.megahit.out"),
            contigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.contigs.fa"),
            fastg = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.fastg"),
            gfa = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.gfa"),
            only_save_scaftigs = \
                config["params"]["assembly"]["megahit"]["only_save_scaftigs"],
            wrapper_dir = WRAPPER_DIR
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.megahit.log")
        script:
            "../wrappers/megahit_wrapper.py"


    rule assembly_megahit_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.gfa.gz")],
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_megahit_all:
        input:


if "idba_ud" in ASSEMBLERS:
    rule assembly_idba_ud:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.idba_ud.out/{assembly_group}.idba_ud.scaftigs.fa.gz"))
        conda:
            config["envs"]["idbaud"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/idba_ud/{assembly_group}.idba_ud.benchmark.txt")
        priority:
            20
        params:
            prefix = "{assembly_group}",
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.idba_ud.out"),
            mink = config["params"]["assembly"]["idba_ud"]["mink"],
            maxk = config["params"]["assembly"]["idba_ud"]["maxk"],
            step = config["params"]["assembly"]["idba_ud"]["step"],
            min_contig = config["params"]["assembly"]["idba_ud"]["min_contig"],
            only_save_scaftigs = \
                config["params"]["assembly"]["idba_ud"]["only_save_scaftigs"]
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.idba_ud.log")
        script:
            "../wrappers/idbaud_wrapper.py"


    rule assembly_idba_ud_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.idba_ud.out/{assembly_group}.idba_ud.scaftigs.fa.gz"),
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_idba_ud_all:
        input:


if "metaspades" in ASSEMBLERS:
    rule assembly_metaspades:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.gfa.gz"))
        conda:
            config["envs"]["spades"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/metaspades/{assembly_group}.metaspades.benchmark.txt")
        priority:
            20
        params:
            prefix = "{assembly_group}",
            memory = str(config["params"]["assembly"]["metaspades"]["memory"]),
            kmers = "auto" \
                if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.metaspades.out"),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["metaspades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = \
                config["params"]["assembly"]["metaspades"]["only_save_scaftigs"],
            link_scaffolds = \
                config["params"]["assembly"]["metaspades"]["link_scaffolds"],
            tar_results = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.tar")
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.metaspades.log")
        script:
            "../wrappers/metaspades_wrapper.py"


    rule assembly_metaspades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.gfa.gz")],
                    assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_metaspades_all:
        input:


if "spades" in ASSEMBLERS:
    rule assembly_spades:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.gfa.gz"))
        conda:
            config["envs"]["spades"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/spades/{assembly_group}.spades.benchmark.txt")
        priority:
            20
        params:
            prefix = "{assembly_group}",
            memory = config["params"]["assembly"]["spades"]["memory"],
            kmers = "auto" \
                if len(config["params"]["assembly"]["spades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["spades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.spades.out"),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["spades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = config["params"]["assembly"]["spades"]["only_save_scaftigs"],
            link_scaffolds = config["params"]["assembly"]["spades"]["link_scaffolds"],
            tar_results = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.tar")
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"], "logs/{assembly_group}.spades.log")
        script:
            "../wrappers/spades_wrapper.py"


    rule assembly_spades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.gfa.gz")],
                    assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_spades_all:
        input:


if "plass" in ASSEMBLERS:
    rule assembly_plass:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            proteins = os.path.join(
                config["output"]["assembly"],
                "proteins/{assembly_group}.plass.out/{assembly_group}.plass.proteins.fa.gz"),
            tmp = directory(temp(os.path.join(
                config["output"]["assembly"],
                "proteins/{assembly_group}.plass.out.tmp")))
        conda:
            config["envs"]["plass"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/plass/{assembly_group}.plass.benchmark.txt")
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.plass.log")
        threads:
            config["params"]["assembly"]["threads"]
        params:
            min_seq_id = config["params"]["assembly"]["plass"]["min_seq_id"],
            min_length = config["params"]["assembly"]["plass"]["min_length"],
            evalue = config["params"]["assembly"]["plass"]["evalue"],
            filter_proteins = config["params"]["assembly"]["plass"]["filter_proteins"]
        shell:
            '''
            plass assemble \
            {input} \
            --threads {threads} \
            --compressed 1 \
            --min-seq-id {params.min_seq_id} \
            --min-length {params.min_length} \
            -e {params.evalue} \
            --filter-proteins {params.filter_proteins} \
            {output.tmp}
            '''


    rule assembly_plass_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "proteins/{assembly_group}.plass.out/{assembly_group}.plass.proteins.fa.gz"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_plass_all:
        input:


def opera_ms_scaftigs_input(wildcards):
    return expand(
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.megahit.scaftigs.fa.gz"),
        sample=wildcards.sample,
        assembler=config["params"]["assembly"]["opera_ms"]["short_read_assembler"])


if "opera_ms" in ASSEMBLERS:
    rule assembly_opera_ms:
        input:
            reads = assembly_input_with_short_and_long_reads,
            scaftigs = opera_ms_scaftigs_input
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.opera_ms.out/{assembly_group}.opera_ms.scaftigs.fa.gz"))
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/opera_ms/{assembly_group}.opera_ms.benchmark.txt")
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.opera_ms.log")
        params:
            opera_ms = config["params"]["assembly"]["opera_ms"]["path"],
            prefix = "{assembly_group}",
            out_dir = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{assembly_group}.opera_ms.out"),
            no_ref_clustering = "--no-ref-clustering" \
                if config["params"]["assembly"]["opera_ms"]["no_ref_clustering"] else "",
            no_strain_clustering = "--no-strain-clustering" \
                if config["params"]["assembly"]["opera_ms"]["no_strain_clustering"] else "",
            no_gap_filling = "--no-gap-filling" \
                if config["params"]["assembly"]["opera_ms"]["no_gap_filling"] else "",
            polishing = "--polishing" \
                if config["params"]["assembly"]["opera_ms"]["polishing"] else "",
            long_read_mapper = config["params"]["assembly"]["opera_ms"]["long_read_mapper"],
            short_read_assembler = config["params"]["assembly"]["opera_ms"]["short_read_assembler"],
            contig_len_threshold = config["params"]["assembly"]["opera_ms"]["contig_len_threshold"],
            contig_edge_len = config["params"]["assembly"]["opera_ms"]["contig_edge_len"],
            contig_window_len = config["params"]["assembly"]["opera_ms"]["contig_window_len"],
            genome_db = "--genome-db %s" % config["params"]["assembly"]["opera_ms"]["genome_db"] \
                if not config["params"]["assembly"]["opera_ms"]["no_ref_clustering"] \
                   else ""
        threads:
            config["params"]["assembly"]["threads"]
        run:
            shell(
                '''
                rm -rf {params.out_dir}
                pigz -p {threads} -dc {input.scaftigs} > {input.scaftigs}.fa

                perl {params.opera_ms} \
                --short-read1 {input.reads[0]} \
                --short-read2 {input.reads[1]} \
                --long-read {input.reads[2]} \
                --contig-file {input.scaftigs}.fa \
                --num-processors {threads} \
                --out-dir {params.out_dir} \
                {params.no_ref_clustering} \
                {params.no_strain_clustering} \
                {params.no_gap_filling} \
                {params.polishing} \
                {params.genome_db} \
                --long-read-mapper {params.long_read_mapper} \
                --short-read-assembler {params.short_read_assembler} \
                --contig-len-thr {params.contig_len_threshold} \
                --contig-edge-len {params.contig_edge_len} \
                --contig-window-len {params.contig_window_len} \
                >{log} 2>&1
                ''')

            shell('''pigz -p {threads} {params.out_dir}/contigs.fasta''')

            if params.polishing != "":
                shell(
                    '''
                    pigz -p {threads} {params.out_dir}/contigs.polished.fasta

                    pushd {params.out_dir} && \
                    ln -s contigs.polished.fasta.gz {params.prefix}.opera_ms.scaftigs.fa.gz && \
                    popd
                    ''')
            else:
                shell(
                    '''
                    pushd {params.out_dir} && \
                    ln -s contigs.fasta.gz {params.prefix}.opera_ms.scaftigs.fa.gz && \
                    popd
                    ''')
            shell('''rm -rf {input.scaftigs}.fa''')


    rule assembly_opera_ms_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.opera_ms.out/{assembly_group}.opera_ms.scaftigs.fa.gz"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_opera_ms_all:
        input:


if len(ASSEMBLERS) != 0:
    if config["params"]["assembly"]["metaquast"]["do"] and IS_PE:
        rule assembly_metaquast:
            input:
                reads = assembly_input_with_short_reads,
                scaftigs = os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
            output:
                protected(os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{assembly_group}.{assembler}.metaquast.out/combined_reference/report.tsv"))
            conda:
                config["envs"]["quast"]
            benchmark:
                os.path.join(config["output"]["assembly"],
                             "benchmark/metaquast_{assembler}/{assembly_group}.{assembler}.metaquast.benchmark.txt")

            log:
                os.path.join(config["output"]["assembly"],
                             "logs/{assembly_group}.{assembler}.metaquast.log")
            params:
                labels = "{assembly_group}.{assembler}",
                output_dir = os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{assembly_group}.{assembler}.metaquast.out")
            threads:
                config["params"]["assembly"]["metaquast"]["threads"]
            shell:
                '''
                metaquast.py \
                {input.scaftigs} \
                --pe1 {input.reads[0]} \
                --pe2 {input.reads[1]} \
                --output-dir {params.output_dir} \
                --labels {params.labels} \
                --circos \
                --ran-finding \
                --conserved-genes-finding \
                --threads {threads} \
                2> {log}
                '''


        rule assembly_metaquast_multiqc:
            input:
                expand(os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{assembly_group}.{{assembler}}.metaquast.out/combined_reference/report.tsv"),
                    assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
            output:
                html = os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                data_dir = directory(
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report_data"))
            conda:
                config["envs"]["multiqc"]
            log:
                os.path.join(config["output"]["assembly"], "logs/multiqc_{assembler}_metaquast.log")
            params:
                output_dir = os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast")
            shell:
                '''
                multiqc \
                --outdir {params.output_dir} \
                --title metaquast \
                --module quast \
                {input} \
                2> {log}
                '''


        rule assembly_metaquast_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["assembly"],
                        "metaquast/{assembly_group}.{assembler}.metaquast.out/combined_reference/report.tsv"),
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report_data")],
                       assembler=ASSEMBLERS,
                       assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

    else:
        rule assembly_metaquast_all:
            input:

           
    rule assembly_report:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            report = os.path.join(
                config["output"]["assembly"],
                "report/{assembler}_stats/{assembly_group}.{assembler}.scaftigs.seqtk.comp.tsv.gz")
        conda:
            config["envs"]["report"]
        priority:
            25
        params:
            assembly_group = "{assembly_group}",
            assembler = "{assembler}"
        shell:
            '''
            seqtk comp {input.scaftigs} | \
            awk \
            'BEGIN \
            {{print "assembly_group\tassembler\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts"}}; \
            {{print "{params.assembly_group}" "\t" "{params.assembler}" "\t" $0}}' | \
            gzip -c > {output.report}
            '''


    rule assembly_report_merge:
        input:
            comp_list = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "report/{{assembler}}_stats/{assembly_group}.{{assembler}}.scaftigs.seqtk.comp.tsv.gz"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
        output:
            summary = os.path.join(
                config["output"]["assembly"],
                "report/assembly_stats_{assembler}.tsv")
        params:
            min_length = config["params"]["assembly"]["report"]["min_length"],
            len_ranges = config["params"]["assembly"]["report"]["len_ranges"]
        threads:
            config["params"]["assembly"]["threads"]
        run:
            comp_list = [(i, params.min_length) for i in input.comp_list]
            metapi.assembler_init(params.len_ranges, ["assembly_group", "assembler"])
            metapi.merge(comp_list, metapi.parse_assembly,
                         threads, output=output.summary)


    rule assembly_report_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "report/assembly_stats_{assembler}.tsv"),
                   assembler=ASSEMBLERS)

else:
    rule assembly_report_all:
        input:


    rule assembly_metaquast_all:
        input:


rule assembly_all:
    input:
        rules.assembly_megahit_all.input,
        rules.assembly_idba_ud_all.input,
        rules.assembly_metaspades_all.input,
        rules.assembly_spades_all.input,
        rules.assembly_plass_all.input,
        rules.assembly_opera_ms_all.input,

        rules.assembly_report_all.input,
        rules.assembly_metaquast_all.input#,

        #rules.rmhost_all.input,
        #rules.qcreport_all.input


localrules:
    assembly_spades_all,
    assembly_idba_ud_all,
    assembly_megahit_all,
    assembly_metaspades_all,
    assembly_opera_ms_all,
    assembly_plass_all,
    assembly_metaquast_all,
    assembly_report_all,
    assembly_all