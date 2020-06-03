def coassembly_inputs(wildcards):
    if RMHOST_DO:
        return get_reads_list("rmhost")
    elif TRIMMING_DO:
        return get_reads_list("trimming")
    else:
        return get_reads_list("raw")


if config["params"]["coassembly"]["megahit"]["do"]:
    rule coassembly_megahit:
        input:
            unpack(coassembly_inputs)
        output:
            scaftigs = os.path.join(config["output"]["coassembly"],
                                    "scaftigs/final.contigs.fa.gz")
        log:
            os.path.join(config["output"]["coassembly"],
                         "logs/coassembly.megahit.log")
        params:
            min_contig = config["params"]["coassembly"]["megahit"]["min_contig"],
            only_save_scaftigs = \
                config["params"]["coassembly"]["megahit"]["only_save_scaftigs"],
            output_dir = os.path.join(config["output"]["coassembly"], "scaftigs")
        threads:
            config["params"]["coassembly"]["megahit"]["threads"]
        run:
            shell('''rm -rf {params.output_dir}''')
            reads_num = len(input)
            if IS_PE:
                shell(
                    '''
                    megahit \
                    -1 %s \
                    -2 %s \
                    -t {threads} \
                    --min-contig-len {params.min_contig} \
                    --out-dir {params.output_dir} \
                    2> {log}
                    ''' % (",".join(input[0:reads_num//2]),
                           ",".join(input[reads_num//2:])))
            else:
                shell(
                    '''
                    megahit \
                    -r %s \
                    -t {threads} \
                    --min-contig-len {params.min_contig} \
                    --out-dir {params.output_dir} \
                    2> {log}
                    ''' % ",".join(input))
            shell('''pigz -p {threads} {params.output_dir}/final.contigs.fa''')

            if params.only_save_scaftigs:
                shell(
                    '''
                    find {params.output_dir} \
                    -type f \
                    ! -wholename "{output.scaftigs}" -delete
                    ''')
                shell('''rm -rf {params.output_dir}/intermediate_contigs''')


    rule coassembly_megahit_all:
        input:
            os.path.join(config["output"]["coassembly"],
                         "scaftigs/final.contigs.fa.gz")

else:
    rule coassembly_megahit_all:
        input:


rule coassembly_all:
    input:
        rules.coassembly_megahit_all.input,

        rules.rmhost_all.input,
        rules.qcreport_all.input
