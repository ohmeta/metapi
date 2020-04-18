def assembly_input(wildcards):
    if RMHOST_DO:
        return get_reads(wildcards, "rmhost", False)
    elif TRIMMING_DO:
        return get_reads(wildcards, "trimming", False)
    else:
        return get_reads(wildcards, "raw", False)

   
if "megahit" in ASSEMBLERS:
    rule assembly_megahit:
        input:
            reads = assembly_input
        output:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.megahit.scaftigs.fa.gz")
        params:
            min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
            output_prefix = "{sample}",
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.megahit.out"),
            contigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.megahit.out/{sample}.contigs.fa"),
            only_save_scaftigs = \
                config["params"]["assembly"]["megahit"]["only_save_scaftigs"]
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{sample}.megahit.log")
        run:
            shell("rm -rf {params.output_dir}")
            if IS_PE:
                shell(
                    '''
                    megahit \
                    -1 {input.reads[0]} \
                    -2 {input.reads[1]} \
                    -t {threads} \
                    --min-contig-len {params.min_contig} \
                    --out-dir {params.output_dir} \
                    --out-prefix {params.output_prefix} \
                    2> {log}
                    ''')
            else:
                shell(
                    '''
                    megahit \
                    -r {input.reads[0]} \
                    -t {threads} \
                    --min-contig-len {params.min_contig} \
                    --out-dir {params.output_dir} \
                    --out-prefix {params.output_prefix} \
                    2> {log}
                    ''')

            shell('''sed -i 's#^>#>{params.output_prefix}_#g' {params.contigs}''')
            shell('''gzip {params.contigs}''')
            shell('''mv {params.contigs}.gz {output.scaftigs}''')

            if params.only_save_scaftigs:
                shell(
                    '''
                    find {params.output_dir} \
                    -type f \
                    ! -wholename "{output.scaftigs}" -delete
                    ''')
                shell('''rm -rf {params.output_dir}/intermediate_contigs''')


    rule assembly_megahit_all:
        input:
            scaftigs = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "scafitgs/{sample}.megahit.out/{sample}.megahit.scaftigs.fa.gz"),
                sample=SAMPLES.index.unique())


if "idba_ud" in ASSEMBLERS:
    rule assembly_idba_ud:
        input:
            reads = assembly_input
        output:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.idba_ud.out/{sample}.idba_ud.scaftigs.fa.gz")
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
                shell('''seqtk seq -A {intput.reads[0]} > %s''' % reads)

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
            scaftigs = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "scafitgs/{sample}.idba_ud.out/{sample}.idba_ud.scaftigs.fa.gz"),
                sample=SAMPLES.index.unique())


metaspades_kmer_list = ["21", "33", "55"]
if len(config["params"]["assembly"]["metaspades"]["kmers"]) > 0:
    metaspades_kmer_list = config["params"]["assembly"]["metaspades"]["kmers"]

spades_kmer_list = ["21", "33", "55"]
if len(config["params"]["assembly"]["spades"]["kmers"]) > 0:
    spades_kmer_list = config["params"]["assembly"]["spades"]["kmers"]


def get_kmer_dirs(wildcards, assembler, kmer_list):
    return expand(
        os.path.join(config["output"]["assembly"],
                     "scaftigs/{sample}.{assembler}.out/K{kmer}"),
        assembler=assembler,
        sample=wildcards.sample,
        kmer=kmer_list)


if "metaspades" in ASSEMBLERS:
    rule assembly_metaspades:
        input:
            reads = assembly_input
        output:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.metaspades.out/{sample}.metaspades.scaftigs.fa.gz")
        params:
            prefix = "{sample}",
            kmers = "auto" \
                if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.metaspades.out"),
            corrected = os.path.join(config["output"]["assembly"],
                                     "scaftigs/{sample}.metaspades.out/corrected"),
            kmer_dirs = lambda wildcards: get_kmer_dirs(wildcards,
                                                        "metaspades",
                                                        metaspades_kmer_list),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["metaspades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = \
                config["params"]["assembly"]["metaspades"]["only_save_scaftigs"],
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
                shell(
                    '''
                    metaspades.py \
                    -1 {input.reads[0]} \
                    -2 {input.reads[1]} \
                    -k {params.kmers} \
                    {params.only_assembler} \
                    --threads {threads} \
                    -o {params.output_dir} \
                    > {log}
                    ''')

                shell('''sed -i 's#^>#>{params.prefix}_#g' {params.output_dir}/scaffolds.fasta''')
                shell('''pigz -p {threads} {params.output_dir}/scaffolds.fasta''')
                shell('''mv {params.output_dir}/scaffolds.fasta.gz {output.scaftigs}''')
                shell('''rm -rf {params.kmer_dirs}''')
                shell('''rm -rf {params.corrected}''')

                if params.only_save_scaftigs:
                    shell(
                        '''
                        find {params.output_dir} \
                        -type f \
                        ! -wholename "{output.scaftigs}" -delete
                        ''')
                else:
                    shell(
                        '''
                        find {params.output_dir} \
                        -type f \
                        ! -wholename "{output.scaftigs}" \
                        ! -wholename "{params.tar_results}" | \
                        xargs -I % sh -c 'tar -rf {params.tar_results} %; rm -rf %'
                        ''')
                    shell('''pigz -p {threads} {params.tar_results}''')

                shell("rm -rf {params.output_dir}/tmp")
                shell("rm -rf {params.output_dir}/misc")
            else:
                print(
                    '''
                    Don't support single-end reads assembly using metaspades\n,
                    you can try SPAdes or MegaHit, IDBA_UD
                    ''')
                sys.exit(1)


    rule assembly_metaspades_all:
        input:
            scaftigs = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "scafitgs/{sample}.metaspades.out/{sample}.metaspades.scaftigs.fa.gz"),
                sample=SAMPLES.index.unique())
 
           
if "spades" in ASSEMBLERS:
    rule assembly_spades:
        input:
            reads = assembly_input
        output:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.spades.out/{sample}.spades.scaftigs.fa.gz")
        params:
            prefix = "{sample}",
            kmers = "auto" \
                if len(config["params"]["assembly"]["spades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["spades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{sample}.spades.out"),
            corrected = os.path.join(config["output"]["assembly"],
                                     "scaftigs/{sample}.spades.out/corrected"),
            kmer_dirs = lambda wildcards: get_kmer_dirs(wildcards, "spades", spades_kmer_list),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["spades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = config["params"]["assembly"]["spades"]["only_save_scaftigs"],
            tar_results = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{sample}.spades.out/{sample}.spades.tar")
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"], "logs/{sample}.spades.log")
        run:
            if IS_PE:
                shell(
                    '''
                    spades.py \
                    -1 {input.reads[0]} \
                    -2 {input.reads[1]} \
                    -k {params.kmers} \
                    {params.only_assembler} \
                    --threads {threads} \
                    -o {params.output_dir} \
                    > {log}
                    ''')
            else:
                shell(
                    '''
                    spades.py \
                    -s {input.reads[0]} \
                    -k {params.kmers} \
                    {params.only_assembler} \
                    --threads {threads} \
                    -o {params.output_dir} \
                    > {log}
                    ''')

            shell('''sed -i 's#^>#>{params.prefix}_#g' {params.output_dir}/scaffolds.fasta''')
            shell('''pigz -p {threads} {params.output_dir}/scaffolds.fasta''')
            shell('''mv {params.output_dir}/scaffolds.fasta.gz {output.scaftigs}''')
            shell('''rm -rf {params.kmer_dirs}''')
            shell('''rm -rf {params.corrected}''')

            if params.only_save_scaftigs:
                shell(
                    '''
                    find {params.output_dir} \
                    -type f \
                    ! -wholename "{output.scaftigs}" -delete
                    ''')
            else:
                shell(
                    '''
                    find {params.output_dir} \
                    -type f \
                    ! -wholename "{output.scaftigs}" \
                    ! -wholename "{params.tar_results}" | \
                    xargs -I % sh -c 'tar -rf {params.tar_results} %; rm -rf %'
                    ''')
                shell('''pigz -p {threads} {params.tar_results}''')

            shell('''rm -rf {params.output_dir}/tmp''')
            shell('''rm -rf {params.output_dir}/misc''')


    rule assembly_spades_all:
        input:
            scaftigs = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "scafitgs/{sample}.spades.out/{sample}.spades.scaftigs.fa.gz"),
                sample=SAMPLES.index.unique())


if len(ASSEMBLERS) != 0:
    rule assembly_report:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            report = os.path.join(
                config["output"]["assembly"],
                "report/{assembler}_stats/{sample}.{assembler}.scaftigs.seqtk.comp.tsv.gz")
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
            metapi.assembler_init(params.len_ranges)
            metapi.merge(input.comp_list, metapi.parse_assembly,
                         threads, save=True, output=output.summary)


    rule assembly_report_all:
        input:
            expand(
                os.path.join(config["output"]["assembly"],
                             "report/assembly_stats_{assembler}.tsv"),
                assembler=ASSEMBLERS
            )


    rule assembly_all:
        input:
            summaries = rules.assembly_report_all.input,
            scaftigs = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz"),
                assembler=ASSEMBLERS,
                sample=SAMPLES.index.unique()
            )

else:
    rule assembly_all:
        input:
