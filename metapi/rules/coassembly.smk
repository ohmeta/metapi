def coassembly_inputs_with_short_reads(wildcards):
    if RMHOST_DO:
        return get_short_reads_list("rmhost")
    elif TRIMMING_DO:
        return get_short_reads_list("trimming")
    else:
        return get_short_reads_list("raw")


if config["params"]["coassembly"]["megahit"]["do"]:
    rule coassembly_megahit:
        input:
            unpack(coassembly_inputs_with_short_reads)
        output:
            scaftigs = os.path.join(
                config["output"]["coassembly"],
                "scaftigs/all.megahit.out/all.megahit.scaftigs.fa.gz"),
            gfa = os.path.join(
                config["output"]["coassembly"],
                "scaftigs/all.megahit.out/all.megahit.scaftigs.gfa.gz")
        priority:
            20
        log:
            os.path.join(config["output"]["coassembly"],
                         "logs/all.megahit.log")
        params:
            output_prefix = "all",
            min_contig = config["params"]["coassembly"]["megahit"]["min_contig"],
            k_list = ",".join(config["params"]["coassembly"]["megahit"]["k_list"]),
            presets = config["params"]["coassembly"]["megahit"]["presets"],
            output_dir = os.path.join(config["output"]["coassembly"],
                                      "scaftigs/all.megahit.out"),
            contigs = os.path.join(
                config["output"]["coassembly"],
                "scaftigs/all.megahit.out/all.contigs.fa"),
            fastg = os.path.join(
                config["output"]["coassembly"],
                "scaftigs/all.megahit.out/all.megahit.scaftigs.fastg"),
            gfa = os.path.join(
                config["output"]["coassembly"],
                "scaftigs/all.megahit.out/all.megahit.scaftigs.gfa"),
            only_save_scaftigs = \
                config["params"]["coassembly"]["megahit"]["only_save_scaftigs"]
        threads:
            config["params"]["coassembly"]["megahit"]["threads"]
        run:
            from Bio import SeqIO
            import re

            if os.path.exists(os.path.join(params.output_dir, "options.json")):
                shell('''megahit --continue --out-dir {params.output_dir}''')
            else:
                reads_num = len(input)
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
                    ''' % ("-1 %s -2 %s" % (",".join(input[0:reads_num//2]),
                                            ",".join(input[reads_num//2:])) \
                           if IS_PE else "-r %s" % ",".join(input),
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


    rule coassembly_megahit_all:
        input:
            os.path.join(config["output"]["coassembly"],
                         "scaftigs/all.megahit.out/all.megahit.scaftigs.fa.gz"),
            os.path.join(config["output"]["coassembly"],
                         "scaftigs/all.megahit.out/all.megahit.scaftigs.gfa.gz")

else:
    rule coassembly_megahit_all:
        input:


rule coassembly_all:
    input:
        rules.coassembly_megahit_all.input,

        rules.rmhost_all.input,
        rules.qcreport_all.input


rule assembly_all:
    input:
        rules.single_assembly_all.input,
        rules.coassembly_all.input
