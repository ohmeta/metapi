def assembly_input_with_short_reads(wildcards):
    if RMHOST_DO:
        return get_reads(wildcards, "rmhost", False)
    elif TRIMMING_DO:
        return get_reads(wildcards, "trimming", False)
    else:
        return get_reads(wildcards, "raw", False)


def assembly_input_with_short_and_long_reads(wildcards):
    if RMHOST_DO:
        return get_reads(wildcards, "rmhost", False, True)
    elif TRIMMING_DO:
        return get_reads(wildcards, "trimming", False, True)
    else:
        return get_reads(wildcards, "raw", False, True)

   
if "megahit" in ASSEMBLERS:
    rule assembly_megahit:
        input:
            reads = assembly_input_with_short_reads
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.gfa.gz"
            ))
        priority:
            20
        params:
            output_prefix = "{sample}",
            min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
            k_list = ",".join(config["params"]["assembly"]["megahit"]["k_list"]),
            presets = config["params"]["assembly"]["megahit"]["presets"],
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.megahit.out"),
            contigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.contigs.fa"),
            fastg = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.fastg"),
            gfa = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.gfa"),
            only_save_scaftigs = \
                config["params"]["assembly"]["megahit"]["only_save_scaftigs"]
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{sample}.megahit.log")
        run:
            from Bio import SeqIO
            import re

            if os.path.exists(os.path.join(params.output_dir, "options.json")):
                shell('''megahit --continue --out-dir {params.output_dir}''')
            else:
                shell("rm -rf {params.output_dir}")
                shell(
                    '''
                    megahit \
                    %s \
                    -t {threads} \
                    %s \
                    --min-contig-len {params.min_contig} \
                    --out-dir {params.output_dir} \
                    --out-prefix {params.output_prefix} \
                    2> {log}
                    ''' % ("-1 {input.reads[0]} -2 {input.reads[1]}" if IS_PE \
                           else "-r {input.reads[0]}",
                           "--presets %s" % params.presets \
                           if params.presets != "" \
                           else "--k-list {params.k_list}"))

            k_num = 0
            for seq_record in SeqIO.parse(params.contigs, "fasta"):
                k_num = int(re.search('k(.*)_', seq_record.id).group(1))
                break

            shell(
                '''
                megahit_toolkit contig2fastg \
                %d \
                {params.contigs} \
                > {params.fastg}
                ''' % k_num)

            shell('''fastg2gfa {params.fastg} > {params.gfa}''')
            shell('''pigz -p {threads} {params.fastg}''')
            shell('''pigz -p {threads} {params.gfa}''')

            shell('''pigz -p {threads} {params.contigs}''')
            shell('''mv {params.contigs}.gz {output.scaftigs}''')

            if params.only_save_scaftigs:
                shell('''fd -t f -E "*.gz" . {params.output_dir} -x rm -rf {{}}''')
                shell('''rm -rf {params.output_dir}/intermediate_contigs''')


    rule assembly_megahit_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.gfa.gz")],
                   sample=SAMPLES.index.unique())

else:
    rule assembly_megahit_all:
        input:


if "idba_ud" in ASSEMBLERS:
    rule assembly_idba_ud:
        input:
            reads = assembly_input_with_short_reads
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.idba_ud.out/{sample}.idba_ud.scaftigs.fa.gz"))
        priority:
            20
        params:
            prefix = "{sample}",
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.idba_ud.out"),
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
                         "logs/{sample}.idba_ud.log")
        run:
            shell('''rm -rf {params.output_dir}''')
            shell('''mkdir {params.output_dir}''')

            reads = os.path.join(
                config["output"]["assembly"],
                "scaftigs/%s.idba_ud.out/%s.fa" % (params.prefix, params.prefix))

            if IS_PE:
                shell(
                    '''
                    seqtk mergepe {input.reads[0]} {input.reads[1]} | \
                    seqtk seq -A - > %s
                    ''' % reads)
            else:
                shell('''seqtk seq -A {input.reads[0]} > %s''' % reads)

            shell(
                '''
                idba_ud \
                -r %s \
                --mink {params.mink} \
                --maxk {params.maxk} \
                --step {params.step} \
                --min_contig {params.min_contig} \
                -o {params.output_dir} \
                --num_threads {threads} \
                --pre_correction \
                > {log}
                ''' % reads)

            shell('''rm -rf %s''' % reads)
            shell('''sed -i 's#^>#>{params.prefix}_#g' {params.output_dir}/scaffold.fa''')
            shell('''pigz -p {threads} {params.output_dir}/scaffold.fa''')
            shell('''mv {params.output_dir}/scaffold.fa.gz {output.scaftigs}''')

            if params.only_save_scaftigs:
                shell(
                    '''
                    find {params.output_dir} \
                    -type f \
                    ! -wholename "{output.scaftigs}" -delete
                    ''')


    rule assembly_idba_ud_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.idba_ud.out/{sample}.idba_ud.scaftigs.fa.gz"),
                   sample=SAMPLES.index.unique())

else:
    rule assembly_idba_ud_all:
        input:


if "metaspades" in ASSEMBLERS:
    rule assembly_metaspades:
        input:
            reads = assembly_input_with_short_reads
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.metaspades.out/{sample}.metaspades.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.metaspades.out/{sample}.metaspades.scaftigs.gfa.gz"))
        priority:
            20
        params:
            prefix = "{sample}",
            memory = config["params"]["assembly"]["metaspades"]["memory"],
            kmers = "auto" \
                if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.metaspades.out"),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["metaspades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = \
                config["params"]["assembly"]["metaspades"]["only_save_scaftigs"],
            link_scaffolds = \
                config["params"]["assembly"]["metaspades"]["link_scaffolds"],
            tar_results = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.metaspades.out/{sample}.metaspades.tar")
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{sample}.metaspades.log")
        run:
            if IS_PE:
                if os.path.exists(os.path.join(params.output_dir, "params.txt")):
                    shell(
                        '''
                        metaspades.py \
                        --continue \
                        -o {params.output_dir} \
                        > {log}
                        ''')
                else:
                    shell('''rm -rf {params.output_dir}''')
                    shell(
                        '''
                        metaspades.py \
                        -1 {input.reads[0]} \
                        -2 {input.reads[1]} \
                        -k {params.kmers} \
                        {params.only_assembler} \
                        --memory {params.memory} \
                        --threads {threads} \
                        --checkpoints last \
                        -o {params.output_dir} \
                        > {log}
                        ''')

                shell(
                    '''
                    pigz -p {threads} {params.output_dir}/scaffolds.fasta
                    mv {params.output_dir}/scaffolds.fasta.gz \
                    {params.output_dir}/{params.prefix}.metaspades.scaffolds.fa.gz
                    ''')
                shell(
                    '''
                    pigz -p {threads} {params.output_dir}/contigs.fasta
                    mv {params.output_dir}/contigs.fasta.gz \
                    {params.output_dir}/{params.prefix}.metaspades.contigs.fa.gz
                    ''')
                shell(
                    '''
                    pigz -p {threads} {params.output_dir}/contigs.paths
                    mv {params.output_dir}/contigs.paths.gz \
                    {params.output_dir}/{params.prefix}.metaspades.contigs.paths.gz
                    ''')
                shell(
                    '''
                    pigz -p {threads} {params.output_dir}/scaffolds.paths
                    mv {params.output_dir}/scaffolds.paths.gz \
                    {params.output_dir}/{params.prefix}.metaspades.scaffolds.paths.gz
                    ''')
                shell(
                    '''
                    pigz -p {threads} {params.output_dir}/assembly_graph_with_scaffolds.gfa
                    mv {params.output_dir}/assembly_graph_with_scaffolds.gfa.gz {output.gfa}
                    ''')

                if params.link_scaffolds:
                    shell(
                        '''
                        pushd {params.output_dir} && \
                        ln -s {params.prefix}.metaspades.scaffolds.fa.gz \
                        {params.prefix}.metaspades.scaftigs.fa.gz && \
                        ln -s {params.prefix}.metaspades.scaffolds.paths.gz \
                        {params.prefix}.metaspades.scaftigs.paths.gz && \
                        popd
                        ''')
                else:
                     shell(
                        '''
                        pushd {params.output_dir} && \
                        ln -s {params.prefix}.metaspades.contigs.fa.gz \
                        {params.prefix}.metaspades.scaftigs.fa.gz && \
                        ln -s {params.prefix}.metaspades.contigs.paths.gz \
                        {params.prefix}.metaspades.scaftigs.paths.gz && \
                        popd
                        ''')

                if params.only_save_scaftigs:
                    shell('''fd -d 1 -E "*.gz" . {params.output_dir} -x rm -rf {{}}''')
                else:
                    shell(
                        '''
                        rm -rf {params.output_dir}/{corrected, misc, pipeline_state, tmp}
                        tar -czvf {params.tar_results}.gz {params.output_dir}/K*
                        ''')
            else:
                print(
                    '''
                    Don't support single-end reads assembly using MetaSPAdes\n,
                    you can try SPAdes or MegaHit, IDBA_UD
                    ''')
                sys.exit(1)


    rule assembly_metaspades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.metaspades.out/{sample}.metaspades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.metaspades.out/{sample}.metaspades.scaftigs.gfa.gz")],
                   sample=SAMPLES.index.unique())

else:
    rule assembly_metaspades_all:
        input:


if "spades" in ASSEMBLERS:
    rule assembly_spades:
        input:
            reads = assembly_input_with_short_reads
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.spades.out/{sample}.spades.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.spades.out/{sample}.spades.scaftigs.gfa.gz"))
        priority:
            20
        params:
            prefix = "{sample}",
            memory = config["params"]["assembly"]["spades"]["memory"],
            kmers = "auto" \
                if len(config["params"]["assembly"]["spades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["spades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.spades.out"),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["spades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = config["params"]["assembly"]["spades"]["only_save_scaftigs"],
            link_scaffolds = config["params"]["assembly"]["spades"]["link_scaffolds"],
            tar_results = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{sample}.spades.out/{sample}.spades.tar")
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"], "logs/{sample}.spades.log")
        run:
            if os.path.exists(os.path.join(params.output_dir, "params.txt")):
                shell(
                    '''
                    spades.py \
                    --continue \
                    -o {params.output_dir} \
                    > {log}
                    ''')
            else:
                shell('''rm -rf {params.output_dir}''')
                shell(
                    '''
                    spades.py \
                    %s \
                    -k {params.kmers} \
                    {params.only_assembler} \
                    --memory {params.memory} \
                    --threads {threads} \
                    --checkpoints last \
                    -o {params.output_dir} \
                    > {log}
                    ''' % "-1 {input.reads[0]} -2 {input.reads[1]}" if IS_PE \
                    else "-s {input.reads[0]}")

            shell(
                '''
                pigz -p {threads} {params.output_dir}/scaffolds.fasta
                mv {params.output_dir}/scaffolds.fasta.gz \
                {params.output_dir}/{params.prefix}.spades.scaffolds.fa.gz
                ''')
            shell(
                '''
                pigz -p {threads} {params.output_dir}/contigs.fasta
                mv {params.output_dir}/contigs.fasta.gz \
                {params.output_dir}/{params.prefix}.spades.contigs.fa.gz
                ''')
            shell(
                '''
                pigz -p {threads} {params.output_dir}/contigs.paths
                mv {params.output_dir}/contigs.paths.gz \
                {params.output_dir}/{params.prefix}.spades.contigs.paths.gz
                ''')
            shell(
                '''
                pigz -p {threads} {params.output_dir}/scaffolds.paths
                mv {params.output_dir}/scaffolds.paths.gz \
                {params.output_dir}/{params.prefix}.spades.scaffolds.paths.gz
                ''')
            shell(
                '''
                pigz -p {threads} {params.output_dir}/assembly_graph_with_scaffolds.gfa
                mv {params.output_dir}/assembly_graph_with_scaffolds.gfa.gz {output.gfa}
                ''')

            if params.link_scaffolds:
                shell(
                    '''
                    pushd {params.output_dir} && \
                    ln -s {params.prefix}.spades.scaffolds.fa.gz \
                    {params.prefix}.spades.scaftigs.fa.gz && \
                    ln -s {params.prefix}.spades.scaffolds.paths.gz \
                    {params.prefix}.spades.scaftigs.paths.gz && \
                    popd
                    ''')
            else:
                shell(
                    '''
                    pushd {params.output_dir} && \
                    ln -s {params.prefix}.spades.contigs.fa.gz \
                    {params.prefix}.spades.scaftigs.fa.gz && \
                    ln -s {params.prefix}.spades.contigs.paths.gz \
                    {params.prefix}.spades.scaftigs.paths.gz && \
                    popd
                    ''')

            if params.only_save_scaftigs:
                shell('''fd -d 1 -E "*.gz" . {params.output_dir} -x rm -rf {{}}''')
            else:
                shell(
                    '''
                    rm -rf {params.output_dir}/{corrected, misc, pipeline_state, tmp}
                    tar -czvf {params.tar_results}.gz {params.output_dir}/K*
                    ''')


    rule assembly_spades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.spades.out/{sample}.spades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.spades.out/{sample}.spades.scaftigs.gfa.gz")],
                   sample=SAMPLES.index.unique())

else:
    rule assembly_spades_all:
        input:


if "plass" in ASSEMBLERS:
    rule assembly_plass:
        input:
            reads = assembly_input_with_short_reads
        output:
            proteins = os.path.join(
                config["output"]["assembly"],
                "proteins/{sample}.plass.out/{sample}.plass.proteins.fa.gz"),
            tmp = directory(temp(os.path.join(
                config["output"]["assembly"],
                "proteins/{sample}.plass.out.tmp")))
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{sample}.plass.log")
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


    rule assembly_nucleotides_plass_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "proteins/{sample}.plass.out/{sample}.plass.proteins.fa.gz"),
                   sample=SAMPLES.index.unique())

else:
    rule assembly_plass_all:
        input:


def opera_ms_scaftigs_input(wildcards):
    return expand(
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{sample}.{assembler}.out/{sample}.megahit.scaftigs.fa.gz"),
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
                "scaftigs/{sample}.opera_ms.out/{sample}.opera_ms.scaftigs.fa.gz"))
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{sample}.opera_ms.log")
        params:
            opera_ms = config["params"]["assembly"]["opera_ms"]["path"],
            prefix = "{sample}",
            out_dir = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{sample}.opera_ms.out"),
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
                "scaftigs/{sample}.opera_ms.out/{sample}.opera_ms.scaftigs.fa.gz"),
                   sample=SAMPLES.index.unique())

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
                    "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
            output:
                protected(os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{sample}.{assembler}.metaquast.out/combined_reference/report.tsv"))
            log:
                os.path.join(config["output"]["assembly"],
                             "logs/{sample}.{assembler}.metaquast.log")
            params:
                labels = "{sample}.{assembler}",
                output_dir = os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{sample}.{assembler}.metaquast.out")
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
                    "metaquast/{sample}.{{assembler}}.metaquast.out/combined_reference/report.tsv"),
                       sample=SAMPLES.index.unique())
            output:
                html = os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                data_dir = directory(
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report_data"))
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
                        "metaquast/{sample}.{assembler}.metaquast.out/combined_reference/report.tsv"),
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report_data")],
                       assembler=ASSEMBLERS,
                       sample=SAMPLES.index.unique())

    else:
        rule assembly_metaquast_all:
            input:

           
    rule assembly_report:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            report = os.path.join(
                config["output"]["assembly"],
                "report/{assembler}_stats/{sample}.{assembler}.scaftigs.seqtk.comp.tsv.gz")
        priority:
            25
        params:
            sample_id = "{sample}",
            assembler = "{assembler}"
        shell:
            '''
            seqtk comp {input.scaftigs} | \
            awk \
            'BEGIN \
            {{print "sample_id\tassembler\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts"}}; \
            {{print "{params.sample_id}" "\t" "{params.assembler}" "\t" $0}}' | \
            gzip -c > {output.report}
            '''


    rule assembly_report_merge:
        input:
            comp_list = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "report/{{assembler}}_stats/{sample}.{{assembler}}.scaftigs.seqtk.comp.tsv.gz"),
                sample=SAMPLES.index.unique())
        output:
            summary = os.path.join(
                config["output"]["assembly"],
                "report/assembly_stats_{assembler}.tsv")
        params:
            len_ranges = config["params"]["assembly"]["report"]["len_ranges"]
        threads:
            config["params"]["assembly"]["threads"]
        run:
            metapi.assembler_init(params.len_ranges, ["sample_id", "assembler"])
            metapi.merge(input.comp_list, metapi.parse_assembly,
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


rule single_assembly_all:
    input:
        rules.assembly_megahit_all.input,
        rules.assembly_idba_ud_all.input,
        rules.assembly_metaspades_all.input,
        rules.assembly_spades_all.input,
        rules.assembly_plass_all.input,
        rules.assembly_opera_ms_all.input,

        rules.assembly_report_all.input,
        rules.assembly_metaquast_all.input,

        rules.rmhost_all.input,
        rules.qcreport_all.input
