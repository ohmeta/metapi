ASSEMBLY_GROUP = SAMPLES.reset_index().loc[:, ["assembly_group", "binning_group"]].drop_duplicates()

assembly_df_list = []
for assembler in ASSEMBLERS:
    assembly_df = ASSEMBLY_GROUP.copy()
    assembly_df["assembler"] = assembler
    assembly_df_list.append(assembly_df)
ASSEMBLY_GROUPS = pd.concat(assembly_df_list, axis=0)


rule assembly_megahit:
    input:
        lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR)
    output:
        scaftigs = os.path.join(config["output"]["assembly"], "scaftigs/{binning_group}.{assembly_group}.megahit/{binning_group}.{assembly_group}.megahit.scaftigs.fa.gz"),
        gfa = os.path.join(config["output"]["assembly"], "scaftigs/{binning_group}.{assembly_group}.megahit/{binning_group}.{assembly_group}.megahit.scaftigs.gfa.gz")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_megahit/{binning_group}.{assembly_group}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_megahit/{binning_group}.{assembly_group}.txt")
    params:
        reads = lambda wildcards: metapi.get_samples_for_assembly_megahit(wildcards, SAMPLES, SAMPLESDIR),
        prefix = "{binning_group}.{assembly_group}",
        min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
        k_list = ",".join(config["params"]["assembly"]["megahit"]["k_list"]),
        presets = config["params"]["assembly"]["megahit"]["presets"],
        only_save_scaftigs = "yes" if config["params"]["assembly"]["megahit"]["only_save_scaftigs"] else "no"
    priority:
        20
    threads:
        config["params"]["assembly"]["threads"]
    conda:
        config["envs"]["megahit"]
    shell:
        '''
        set +e

        OUTDIR=$(dirname {output.scaftigs})
        CONTIGS=$OUTDIR/{params.prefix}.contigs.fa
        FASTG=$OUTDIR/{params.prefix}.megahit.scaftigs.fastg.gz

        if [ -e $OUTDIR/options.json ];
        then
            megahit --continue --out-dir $OUTDIR
        else
            rm -rf $OUTDIR
            KMEROPTS=""
            if [ "{params.presets}" != "" ];
            then
                KMEROPTS="--presets {params.presets}"
            else
                KMEROPTS="--k-list {params.k_list}"
            fi

            megahit \
            {params.reads} \
            -t {threads} \
            $KMEROPTS \
            --min-contig-len {params.min_contig} \
            --out-dir $OUTDIR \
            --out-prefix {params.output_prefix} \
            >{log} 2>&1
        fi

        if [ -f $CONTIGS ];
        then
            knum=`grep "^>" {params.contigs} | head -1 | sed 's/>k//g' | awk -F_ '{{print $1}}'`
            megahit_toolkit contig2fastg $knum $CONTIGS | pigz -f -p {threads} > $FASTG
            fastg2gfa $FASTG | pigz -f -p {threads} > {output.gfa}

            pigz -f -p {threads} $CONTIGS
            mv $CONTIGS.gz {output.scaftigs}

            if [ "{params.only_save_scaftigs}" == "yes" ];
            then
                fd -t f -E "*.gz" . $OUTDIR -x rm -rf {{}}
                rm -rf $OUTDIR/intermediate_contigs
            fi
        else
            touch $OUTDIR/FAILED
        fi
        '''


if "megahit" in ASSEMBLERS:
    rule assembly_megahit_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.megahit/{binning_group}.{assembly_group}.megahit.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.megahit/{binning_group}.{assembly_group}.megahit.scaftigs.gfa.gz")],
                zip,
                binning_group=ASSEMBLY_GROUP["binning_group"],
                assembly_group=ASSEMBLY_GROUP["assembly_group"])

else:
    rule assembly_megahit_all:
        input:


rule assembly_idba_ud:
    input:
        lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR)
    output:
        os.path.join(config["output"]["assembly"], "scaftigs/{binning_group}.{assembly_group}.idba_ud/{binning_group}.{assembly_group}.idba_ud.scaftigs.fa.gz")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_idba_ud/{binning_group}.{assembly_group}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_idba_ud/{binning_group}.{assembly_group}.txt")
    params:
        reads_cmd = lambda wildcards: metapi.get_samples_for_assembly_idba_ud(wildcards, SAMPLES, SAMPLESDIR),
        samples_dir = SAMPLESDIR,
        prefix = "{binning_group}.{assembly_group}",
        mink = config["params"]["assembly"]["idba_ud"]["mink"],
        maxk = config["params"]["assembly"]["idba_ud"]["maxk"],
        step = config["params"]["assembly"]["idba_ud"]["step"],
        min_contig = config["params"]["assembly"]["idba_ud"]["min_contig"],
        only_save_scaftigs = "yes" if config["params"]["assembly"]["idba_ud"]["only_save_scaftigs"] else "no"
    priority:
        20
    threads:
        config["params"]["assembly"]["threads"]
    conda:
        config["envs"]["idba_ud"]
    shell:
        '''
        OUTDIR=$(dirname {output})
        rm -rf $OUTDIR
        mkdir -p $OUTDIR

        {params.reads_cmd[0]}

        idba_ud \
        {params.reads_cmd[1]} \
        --mink {params.mink} \
        --maxk {params.maxk} \
        --step {params.step} \
        --min_contig {params.min_contig} \
        -o $OUTDIR \
        --num_threads {threads} \
        --pre_correction \
        >{log} 2>&1

        rm -rf {params.reads_cmd[2]}

        sed -i 's#^>#>{params.prefix}_#g' $OUTDIR/scaffold.fa
        pigz -f -p {threads} $OUTDIR/scaffold.fa
        mv $OUTDIR/scaffold.fa.gz {output}

        if [ "{params.only_save_scaftigs}" == "yes" ];
        then
            find $OUTDIR -type f ! -wholename "{output}" -delete
        fi
        '''


if "idba_ud" in ASSEMBLERS:
    rule assembly_idba_ud_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.idba_ud/{binning_group}.{assembly_group}.idba_ud.scaftigs.fa.gz"),
                zip,
                binning_group=ASSEMBLY_GROUP["binning_group"],
                assembly_group=ASSEMBLY_GROUP["assembly_group"])

else:
    rule assembly_idba_ud_all:
        input:


rule assembly_metaspades:
    input:
        lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR)
    output:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.metaspades/{binning_group}.{assembly_group}.metaspades.scaftigs.fa.gz"),
        gfa = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.metaspades/{binning_group}.{assembly_group}.metaspades.scaftigs.gfa.gz")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_metaspades/{binning_group}.{assembly_group}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_metaspades/{binning_group}.{assembly_group}.txt")
    params:
        opts = lambda wildcards: metapi.parse_assembly_spades_params(wildcards, config["output"]["assembly"], "metaspades"),
        dataset = lambda wildcards: metapi.get_samples_for_assembly_spades(wildcards, SAMPLES, SAMPLESDIR, "metaspades"),
        prefix = "{binning_group}.{assembly_group}",
        memory = str(config["params"]["assembly"]["metaspades"]["memory"]),
        kmers = "auto" \
        if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 \
        else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
        only_assembler = "--only-assembler" if config["params"]["assembly"]["metaspades"]["only_assembler"] else "",
        only_save_scaftigs = config["params"]["assembly"]["metaspades"]["only_save_scaftigs"],
        link_scaffolds = "yes" if config["params"]["assembly"]["metaspades"]["link_scaffolds"] else "no"
    priority:
        20
    threads:
        config["params"]["assembly"]["threads"]
    conda:
        config["envs"]["spades"]
    shell:
        '''
        OUTDIR=$(dirname {output.scaftigs})

        errorcorrect="no"
        if [ "{params.only_assembler}" != "" ];
        then
            errorcorrect="yes"
        fi

        if [ "{params.kmers}" == "{params.opts[0]}" ] && [ "{params.memory}" == "{params.opts[1]}" ] && [ "{threads}" == "{params.opts[2]}" ] && [ $errorcorrect == "{params.opts[3]}" ];
        then
            metaspades.py \
            --continue \
            -o $OUTDIR \
            >{log} 2>&1
        else
            rm -rf $OUTDIR
            mkdir -p $OUTDIR

            metaspades.py \
            --dataset {params.dataset} \
            -k {params.kmers} \
            {params.only_assembler} \
            --memory {params.memory} \
            --threads {threads} \
            --checkpoints last \
            -o $OUTDIR \
            >{log} 2>&1
        fi

        pigz -f -p {threads} $OUTDIR/scaffolds.fasta
        mv $OUTDIR/scaffolds.fasta.gz $OUTDIR/{params.prefix}.metaspades.scaffolds.fa.gz

        pigz -f -p {threads} $OUTDIR/contigs.fasta
        mv $OUTDIR/contigs.fasta.gz $OUTDIR/{params.prefix}.metaspades.contigs.fa.gz

        pigz -f -p {threads} $OUTDIR/contigs.paths
        mv $OUTDIR/contigs.paths.gz $OUTDIR/{params.prefix}.metaspades.contigs.paths.gz

        pigz -f -p {threads} $OUTDIR/scaffolds.paths
        mv $OUTDIR/scaffolds.paths.gz $OUTDIR/{params.prefix}.metaspades.scaffolds.paths.gz

        pigz -f -p {threads} $OUTDIR/assembly_graph_with_scaffolds.gfa
        mv $OUTDIR/assembly_graph_with_scaffolds.gfa.gz {output.gfa}

        if [ "{params.link_scaffolds}" == "yes" ];
        then
            pushd $OUTDIR
            ln -s {params.prefix}.metaspades.scaffolds.fa.gz {params.prefix}.metaspades.scaftigs.fa.gz
            ln -s {params.prefix}.metaspades.scaffolds.paths.gz {params.prefix}.metaspades.scaftigs.paths.gz
            popd
        else
            pushd $OUTDIR
            ln -s {params.prefix}.metaspades.contigs.fa.gz {params.prefix}.metaspades.scaftigs.fa.gz
            ln -s {params.prefix}.metaspades.contigs.paths.gz {params.prefix}.metaspades.scaftigs.paths.gz
            popd
        fi

        if [ "{params.only_save_scaftigs}" == "True" ];
        then
            fd -d 1 -E "*.gz" . $OUTDIR -x rm -rf {{}}
        else
            TARRES=$OUTDIR/{params.prefix}.metaspades.tar.gz
            rm -rf $OUTDIR/{{corrected, misc, pipeline_state, tmp}}
            tar -czvf $TARRES $OUTDIR/K*
        fi
        '''


if "metaspades" in ASSEMBLERS:
    rule assembly_metaspades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.metaspades/{binning_group}.{assembly_group}.metaspades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.metaspades/{binning_group}.{assembly_group}.metaspades.scaftigs.gfa.gz")],
                    zip,
                    binning_group=ASSEMBLY_GROUP["binning_group"],
                    assembly_group=ASSEMBLY_GROUP["assembly_group"])

else:
    rule assembly_metaspades_all:
        input:


rule assembly_spades:
    input:
        lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR)
    output:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.spades/{binning_group}.{assembly_group}.spades.scaftigs.fa.gz"),
        gfa = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.spades/{binning_group}.{assembly_group}.spades.scaftigs.gfa.gz")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_spades/{binning_group}.{assembly_group}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_spades/{binning_group}.{assembly_group}.txt")
    params:
        opts = lambda wildcards: metapi.parse_assembly_spades_params(wildcards, config["output"]["assembly"], "spades"),
        dataset = lambda wildcards: metapi.get_samples_for_assembly_spades(wildcards, SAMPLES, SAMPLESDIR, "spades"),
        prefix = "{binning_group}.{assembly_group}",
        memory = str(config["params"]["assembly"]["metaspades"]["memory"]),
        kmers = "auto" \
        if len(config["params"]["assembly"]["spades"]["kmers"]) == 0 \
        else ",".join(config["params"]["assembly"]["spades"]["kmers"]),
        only_assembler = "--only-assembler" if config["params"]["assembly"]["spades"]["only_assembler"] else "",
        only_save_scaftigs = config["params"]["assembly"]["spades"]["only_save_scaftigs"],
        link_scaffolds = "yes" if config["params"]["assembly"]["spades"]["link_scaffolds"] else "no"
    priority:
        20
    threads:
        config["params"]["assembly"]["threads"]
    conda:
        config["envs"]["spades"]
    shell:
        '''
        OUTDIR=$(dirname {output.scaftigs})

        errorcorrect="no"
        if [ "{params.only_assembler}" != "" ];
        then
            errorcorrect="yes"
        fi

        if [ "{params.kmers}" == "{params.opts[0]}" ] && [ "{params.memory}" == "{params.opts[1]}" ] && [ "{threads}" == "{params.opts[2]}" ] && [ $errorcorrect == "{params.opts[3]}" ];
        then
            spades.py \
            --continue \
            -o $OUTDIR \
            >{log} 2>&1
        else
            rm -rf $OUTDIR
            mkdir -p $OUTDIR

            spades.py \
            --dataset {params.dataset} \
            -k {params.kmers} \
            {params.only_assembler} \
            --memory {params.memory} \
            --threads {threads} \
            --checkpoints last \
            -o $OUTDIR \
            >{log} 2>&1
        fi

        pigz -f -p {threads} $OUTDIR/scaffolds.fasta
        mv $OUTDIR/scaffolds.fasta.gz $OUTDIR/{params.prefix}.spades.scaffolds.fa.gz

        pigz -f -p {threads} $OUTDIR/contigs.fasta
        mv $OUTDIR/contigs.fasta.gz $OUTDIR/{params.prefix}.spades.contigs.fa.gz

        pigz -f -p {threads} $OUTDIR/contigs.paths
        mv $OUTDIR/contigs.paths.gz $OUTDIR/{params.prefix}.spades.contigs.paths.gz

        pigz -f -p {threads} $OUTDIR/scaffolds.paths
        mv $OUTDIR/scaffolds.paths.gz $OUTDIR/{params.prefix}.spades.scaffolds.paths.gz

        pigz -f -p {threads} $OUTDIR/assembly_graph_with_scaffolds.gfa
        mv $OUTDIR/assembly_graph_with_scaffolds.gfa.gz {output.gfa}

        if [ "{params.link_scaffolds}" == "yes" ];
        then
            pushd $OUTDIR
            ln -s {params.prefix}.spades.scaffolds.fa.gz {params.prefix}.spades.scaftigs.fa.gz
            ln -s {params.prefix}.spades.scaffolds.paths.gz {params.prefix}.spades.scaftigs.paths.gz
            popd
        else
            pushd $OUTDIR
            ln -s {params.prefix}.spades.contigs.fa.gz {params.prefix}.spades.scaftigs.fa.gz
            ln -s {params.prefix}.spades.contigs.paths.gz {params.prefix}.spades.scaftigs.paths.gz
            popd
        fi

        if [ "{params.only_save_scaftigs}" == "True" ];
        then
            fd -d 1 -E "*.gz" . $OUTDIR -x rm -rf {{}}
        else
            TARRES=$OUTDIR/{params.prefix}.spades.tar.gz
            rm -rf $OUTDIR/{{corrected, misc, pipeline_state, tmp}}
            tar -czvf $TARRES $OUTDIR/K*
        fi
        '''


if "spades" in ASSEMBLERS:
    rule assembly_spades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.spades/{binning_group}.{assembly_group}.spades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.spades/{binning_group}.{assembly_group}.spades.scaftigs.gfa.gz")],
                    zip,
                    binning_group=ASSEMBLY_GROUP["binning_group"],
                    assembly_group=ASSEMBLY_GROUP["assembly_group"])

else:
    rule assembly_spades_all:
        input:


rule assembly_plass:
    input:
        lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR)
    output:
        os.path.join(config["output"]["assembly"], "proteins/{binning_group}.{assembly_group}.plass/done")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_plass/{binning_group}.{assembly_group}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_plass/{binning_group}.{assembly_group}.txt")
    params:
        reads_cmd = lambda wildcards: metapi.get_samples_for_assembly_plass(wildcards, SAMPLES, SAMPLESDIR),
        min_seq_id = config["params"]["assembly"]["plass"]["min_seq_id"],
        min_length = config["params"]["assembly"]["plass"]["min_length"],
        evalue = config["params"]["assembly"]["plass"]["evalue"],
        filter_proteins = config["params"]["assembly"]["plass"]["filter_proteins"]
    priority:
        20
    threads:
        config["params"]["assembly"]["threads"]
    conda:
        config["envs"]["plass"]
    shell:
        '''
        OUTDIR=$(dirname {output})
        OUTPE=$OUTDIR/pe
        OUTSE=$OUTDIR/se

        rm -rf $OUTDIR $OUTPE $OUTSE
        mkdir -p $OUTDIR

        {params.reads_cmd[0]}

        if [ "{params.reads_cmd[1]}" != " " ];
        then
            plass assemble \
            {params.reads_cmd[1]} \
            $OUTPE \
            $OUTPE.temp \
            --threads {threads} \
            --compressed 1 \
            --min-seq-id {params.min_seq_id} \
            --min-length {params.min_length} \
            -e {params.evalue} \
            --filter-proteins {params.filter_proteins} \
            --remove-tmp-files 1 \
            >{log} 2>&1
        fi

        if [ "{params.reads_cmd[2]}" != " " ];
        then
            plass assemble \
            {params.reads_cmd[2]} \
            $OUTSE \
            $OUTSE.temp \
            --threads {threads} \
            --compressed 1 \
            --min-seq-id {params.min_seq_id} \
            --min-length {params.min_length} \
            -e {params.evalue} \
            --filter-proteins {params.filter_proteins} \
            --remove-tmp-files 1 \
            >>{log} 2>&1
        fi

        rm -rf {params.reads_cmd[3]}
        touch {output}
        '''


if "plass" in ASSEMBLERS:
    rule assembly_plass_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "proteins/{binning_group}.{assembly_group}.plass/done"),
                zip,
                binning_group=ASSEMBLY_GROUP["binning_group"],
                assembly_group=ASSEMBLY_GROUP["assembly_group"])

else:
    rule assembly_plass_all:
        input:


rule assembly_opera_ms:
    input:
        samples = lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR),
        scaftigs = os.path.join(config["output"]["assembly"], "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        os.path.join(config["output"]["assembly"], "scaftigs/{binning_group}.{assembly_group}.opera_ms/{binning_group}.{assembly_group}.opera_ms.scaftigs.fa.gz")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_opera_ms/{binning_group}.{assembly_group}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_opera_ms/{binning_group}.{assembly_group}.txt")
    params:
        opera_ms = config["params"]["assembly"]["opera_ms"]["path"],
        prefix = "{binning_group}.{assembly_group}",
        reads_cmd = lambda wildcards: metapi.get_samples_for_assembly_opera_ms(wildcards, SAMPLES, SAMPLESDIR, config["output"]["assembly"]),
        no_ref_clustering = "--no-ref-clustering" if config["params"]["assembly"]["opera_ms"]["no_ref_clustering"] else "",
        no_strain_clustering = "--no-strain-clustering" if config["params"]["assembly"]["opera_ms"]["no_strain_clustering"] else "",
        no_gap_filling = "--no-gap-filling" if config["params"]["assembly"]["opera_ms"]["no_gap_filling"] else "",
        polishing = "--polishing" if config["params"]["assembly"]["opera_ms"]["polishing"] else "",
        long_read_mapper = config["params"]["assembly"]["opera_ms"]["long_read_mapper"],
        short_read_assembler = config["params"]["assembly"]["opera_ms"]["short_read_assembler"],
        contig_len_threshold = config["params"]["assembly"]["opera_ms"]["contig_len_threshold"],
        contig_edge_len = config["params"]["assembly"]["opera_ms"]["contig_edge_len"],
        contig_window_len = config["params"]["assembly"]["opera_ms"]["contig_window_len"],
        genome_db = "--genome-db %s" % config["params"]["assembly"]["opera_ms"]["genome_db"] \
        if not config["params"]["assembly"]["opera_ms"]["no_ref_clustering"] \
        else ""
    priority:
        20
    threads:
        config["params"]["assembly"]["threads"]
    run:
        shell(
            '''
            OUTDIR=$(dirname {output})
            rm -rf $OUTDIR

            {params.reads_cmd[0]}

            perl {params.opera_ms} \
            {params.reads_cmd[1]} \
            --num-processors {threads} \
            --out-dir $OUTDIR \
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

            rm -rf {params.reads_cmd[2]}

            pigz -f -p {threads} $OUTDIR/contigs.fasta

            if [ "{params.polishing}" != "" ];
            then
                pigz -f -p {threads} $OUTDIR/contigs.polished.fasta
                pushd $OUTDIR && \
                ln -s contigs.polished.fasta.gz {params.prefix}.opera_ms.scaftigs.fa.gz && \
                popd
            else
                pushd $OUTDIR && \
                ln -s contigs.fasta.gz {params.prefix}.opera_ms.scaftigs.fa.gz && \
                popd
            fi
            '''


if "opera_ms" in ASSEMBLERS:
    rule assembly_opera_ms_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.opera_ms/{binning_group}.{assembly_group}.opera_ms.scaftigs.fa.gz"),
                zip,
                binning_group=ASSEMBLY_GROUP["binning_group"],
                assembly_group=ASSEMBLY_GROUP["assembly_group"])

else:
    rule assembly_opera_ms_all:
        input:


rule assembly_metaquast:
    input:
        samples = lambda wildcards: metapi.get_samples_for_assembly_list(wildcards, SAMPLES, SAMPLESDIR),
        scaftigs = os.path.join(config["output"]["assembly"], "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        os.path.join(
            config["output"]["assembly"],
            "report/metaquast/{binning_group}.{assembly_group}.{assembler}/combined_reference/report.tsv")
    log:
        os.path.join(
            config["output"]["assembly"],
            "logs/assembly_metaquast/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["assembly"],
            "benchmark/assembly_metaquast/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        labels = "{binning_group}.{assembly_group}.{assembler}",
        reads_cmd = lambda wildcards: metapi.get_samples_for_metaquast(wildcards, SAMPLES, SAMPLESDIR)
    threads:
        config["params"]["assembly"]["metaquast"]["threads"]
    conda:
        config["envs"]["quast"]
    shell:
        '''
        OUTDIR=$(dirname $(dirname {output}))

        {params.reads_cmd[0]}

        metaquast.py \
        {params.reads_cmd[1]} \
        --output-dir $OUTDIR \
        --labels {params.labels} \
        --circos \
        --ran-finding \
        --conserved-genes-finding \
        --threads {threads} \
        >{log} 2>&1

        rm -rf {params.reads_cmd[2]}
        '''


rule assembly_metaquast_multiqc:
    input:
        expand(os.path.join(
            config["output"]["assembly"],
            "metaquast/{binning_group}.{assembly_group}.{{assembler}}.metaquast/combined_reference/report.tsv"),
            binning_group=ASSEMBLY_GROUP["binning_group"],
            assembly_group=ASSEMBLY_GROUP["assembly_group"])
    output:
        os.path.join(config["output"]["assembly"], "report/multiqc/{assembler}/metaquast_multiqc_report.html")
    log:
        os.path.join(config["output"]["assembly"], "logs/assembly_metaquast_multiqc/{assembler}.log")
    benchmark:
        os.path.join(config["output"]["assembly"], "benchmark/assembly_metaquast_multiqc/{assembler}.txt")
    conda:
        config["envs"]["multiqc"]
    threads:
        1
    shell:
        '''
        OUTDIR=$(dirname {output})
        rm -rf $OUTDIR

        multiqc \
        --outdir $OUTDIR \
        --title metaquast \
        --module quast \
        {input} \
        >{log} 2>&1
        '''


if config["params"]["assembly"]["metaquast"]["do"]:
    rule assembly_metaquast_all:
        input:
            expand(
                os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{binning_group}.{assembly_group}.{assembler}.metaquast/combined_reference/report.tsv"),
                    zip,
                    assembler=ASSEMBLY_GROUPS["assembler"],
                    binning_group=ASSEMBLY_GROUPS["binning_group"],
                    assembly_group=ASSEMBLY_GROUPS["assembly_group"]),
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast/metaquast_multiqc_report_data")],
                    assembler=ASSEMBLERS)

else:
    rule assembly_metaquast_all:
        input:


rule assembly_report:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        os.path.join(
            config["output"]["assembly"],
            "report/{assembler}_stats/{binning_group}.{assembly_group}.{assembler}.scaftigs.seqtk.comp.tsv.gz")
    params:
        binning_group = "{binning_group}",
        assembly_group = "{assembly_group}",
        assembler = "{assembler}"
    priority:
        25
    threads:
        1
    conda:
        config["envs"]["report"]
    shell:
        '''
        seqtk comp {input} | \
        awk \
        'BEGIN \
        {{print "binning_group\tassembly_group\tassembler\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts"}}; \
        {{print "{params.binning_group}\t{params.assembly_group}\t{params.assembler}\t" $0}}' | \
        gzip -c > {output}
        '''


rule assembly_report_merge:
    input:
        comp_list = expand(
            os.path.join(config["output"]["assembly"], "report/{{assembler}}_stats/{binning_group}.{assembly_group}.{{assembler}}.scaftigs.seqtk.comp.tsv.gz"),
            zip,
            binning_group=ASSEMBLY_GROUP["binning_group"],
            assembly_group=ASSEMBLY_GROUP["assembly_group"])
    output:
        summary = os.path.join(config["output"]["assembly"], "report/assembly_stats_{assembler}.tsv.gz")
    params:
        min_length = config["params"]["assembly"]["report"]["min_length"],
        len_ranges = config["params"]["assembly"]["report"]["len_ranges"]
    threads:
        config["params"]["assembly"]["threads"]
    run:
        comp_list = [(i, params.min_length) for i in input.comp_list]
        metapi.assembler_init(params.len_ranges, ["binning_group", "assembly_group", "assembler"])
        metapi.merge(comp_list, metapi.parse_assembly, threads, output=output.summary)


rule assembly_report_all:
    input:
        expand(os.path.join(config["output"]["assembly"], "report/assembly_stats_{assembler}.tsv.gz"),
        assembler=ASSEMBLERS)


rule assembly_all:
    input:
        rules.assembly_megahit_all.input,
        rules.assembly_idba_ud_all.input,
        rules.assembly_metaspades_all.input,
        rules.assembly_spades_all.input,
        rules.assembly_plass_all.input,
        rules.assembly_opera_ms_all.input,
        rules.assembly_report_all.input,
        rules.assembly_metaquast_all.input


localrules:
    assembly_megahit_all,
    assembly_idba_ud_all,
    assembly_metaspades_all,
    assembly_spades_all,
    assembly_plass_all,
    assembly_opera_ms_all,
    assembly_metaquast_all,
    assembly_report_all,
    assembly_all
