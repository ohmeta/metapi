if config["params"]["binning"]["graphbin2"]["do"]:
    rule binning_graphbin2_prepare_assembly:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
            gfa = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.gfa.gz")
        output:
             scaftigs = temp(os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/graphbin2/scaftigs.fa")),
             gfa = temp(os.path.join(
                 config["output"]["binning"],
                 "mags/{binning_group}.{assembly_group}.{assembler}/graphbin2/scaftigs.gfa"))
        conda:
            config["envs"]["report"]
        shell:
            '''
            pigz -fdc {input.scaftigs} > {output.scaftigs}
            pigz -fdc {input.gfa} > {output.gfa}
            '''

  
    rule binning_graphbin2_prepare_binned:
        input:
            mags_dir = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/{binner_graphbin}")
        output:
            binned = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/graphbin2/{binning_group}.{assembly_group}.{assembler}.{binner_graphbin}.graphbin2.csv")
        params:
            assembler = "{assembler}"
        run:
            metapi.get_binning_info(input.mags_dir,
                                    output.binned,
                                    params.assembler)


    localrules:
        binning_graphbin2_prepare_binned


    rule binning_graphbin2:
        input:
            scaftigs = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/graphbin2/scaftigs.fa"),
            gfa = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/graphbin2/scaftigs.gfa"),
            binned = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/graphbin2/{binning_group}.{assembly_group}.{assembler}.{binner_graphbin}.graphbin2.csv")
        output:
            os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/{binner_graphbin}_graphbin2/binning_done")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/{binning_graphbin}/{binning_group}.{assembly_group}.{assembler}.{binner_graphbin}.graphbin2.refine.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/{binner_graphbin}/{binning_group}.{assembly_group}.{assembler}.{binner_graphbin}.benchmark.txt")
        params:
            assembler = "{assembler}",
            mags_dir = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/{binner_graphbin}_graphbin2/"),
            prefix = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/{binner_graphbin}_graphbin2/{binning_group}.{assembly_group}.{assembler}.{binner_graphbin}_graphbin2.bin"),
            paths = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.paths.gz"),
            depth = config["params"]["binning"]["graphbin2"]["depth"],
            threshold = config["params"]["binning"]["graphbin2"]["threshold"]
        threads:
            config["params"]["binning"]["threads"]
        run:
            import pandas as pd
            import os

            shell('''mkdir -p {params.mags_dir}''')

            df = pd.read_csv(input.binned, names=["scaftigs_id", "bin_id"])

            if not df.empty:
                if params.assembler == "metaspades" or params.assembler == "spades":
                    shell(
                        '''
                        pigz -f -p {threads} -dc {params.paths} > {params.mags_dir}/scaftigs.paths

                        graphbin2 \
                        --assembler spades \
                        --contigs {input.scaftigs} \
                        --graph {input.gfa} \
                        --paths {params.mags_dir}/scaftigs.paths \
                        --binned {input.binned} \
                        --nthreads {threads} \
                        --depth {params.depth} \
                        --threshold {params.threshold} \
                        --output {params.mags_dir} \
                        > {log} 2>&1

                        rm -rf {params.mags_dir}/scaftigs.paths
                        ''')
                else:
                    shell(
                        '''
                        graphbin2 \
                        --assembler {params.assembler} \
                        --contigs {input.scaftigs} \
                        --graph {input.gfa} \
                        --binned {input.binned} \
                        --nthreads {threads} \
                        --depth {params.depth} \
                        --threshold {params.threshold} \
                        --output {params.mags_dir} \
                        > {log} 2>&1
                        ''')

                metapi.generate_mags(f"{params.mags_dir}/graphbin2_output.csv",
                                     input.scaftigs, params.prefix)
                shell('''touch {output}''')


    rule binning_graphbin2_all:
        input:
            expand(expand(os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/{{binner_graphbin}}_graphbin2/binning_done"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
                binner_graphbin=BINNERS_GRAPHBIN)

            #rules.alignment_all.input,
            #rules.assembly_all.input

else:
    rule binning_graphbin2_all:
        input:


def get_binning_done(wildcards, binners):
    return expand(os.path.join(
        config["output"]["binning"],
        "mags/{binning_group}.{assembly_group}.{assembler}/{binner}/binning_done"),
        binning_group=wildcards.binning_group,
        assembly_group=wildcards.assembly_group,
        assembler=wildcards.assembler,
        binner=binners)


if config["params"]["binning"]["dastools"]["do"]:
    rule binning_dastools_preprocess:
        input:
            lambda wildcards: get_binning_done(wildcards, BINNERS_DASTOOLS)
        output:
            expand(os.path.join(
                config["output"]["binning"],
                "mags_id/{{binning_group}}.{{assembly_group}}.{{assembler}}/{binner_dastools}_Contigs2Bin.tsv"),
                binner_dastools=BINNERS_DASTOOLS)
        run:
            import glob
            import gzip
            import os
            from Bio import SeqIO
            from natsort import natsorted

            i = -1
            for binning_done in input:
                i += 1
                mags_dir = os.path.dirname(binning_done)
                shell(f'''rm -rf {output[i]}''')
                mags_list = glob.glob(mags_dir + "/*.bin.*.fa.gz")

                if len(mags_list) == 0:
                    shell(f'''touch {output[i]}''')
                else:
                    with open(output[i], 'w') as oh:
                        for bin_fa in natsorted(mags_list):
                            bin_id_list = os.path.basename(bin_fa).split(".")
                            bin_id = bin_id_list[-5] + "." + str(bin_id_list[-3])
                            with gzip.open(bin_fa, 'rt') as ih:
                                for contig in SeqIO.parse(ih, "fasta"):
                                    oh.write(f'''{contig.id}\t{bin_id}\n''')


    rule binning_dastools:
        input:
            contigs2bin = expand(
                os.path.join(
                    config["output"]["binning"],
                    "mags_id/{{binning_group}}.{{assembly_group}}.{{assembler}}/{binner_dastools}_Contigs2Bin.tsv"),
                    binner_dastools=BINNERS_DASTOOLS),
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
            pep = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binning_group}.{assembly_group}.{assembler}.prodigal.faa.gz")
        output:
            os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/dastools/binning_done")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/dastools/{binning_group}.{assembly_group}.{assembler}.dastools.binning.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/dastools/{binning_group}.{assembly_group}.{assembler}.dastools.benchmark.txt")
        priority:
            30
        conda:
            config["envs"]["dastools"]
        params:
            binner = ",".join(BINNERS_DASTOOLS),
            wrapper_dir = WRAPPER_DIR,
            mags_dir = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/dastools"),
            search_engine = config["params"]["binning"]["dastools"]["search_engine"],
            score_threshold = config["params"]["binning"]["dastools"]["score_threshold"],
            duplicate_penalty = config["params"]["binning"]["dastools"]["duplicate_penalty"],
            megabin_penalty = config["params"]["binning"]["dastools"]["megabin_penalty"],
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/dastools/{binning_group}.{assembly_group}.{assembler}.dastools.bin")
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            set +e
            rm -rf {params.mags_dir}
            mkdir -p {params.mags_dir}
            
            contigs2bin=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.contigs2bin})

            FNA={input.scaftigs}
            PEP={input.pep}

            pigz -dkf $FNA
            pigz -dkf $PEP

            DAS_Tool \
            --bins $contigs2bin \
            --labels {params.binner} \
            --contigs ${{FNA%.gz}} \
            --proteins ${{PEP%.gz}} \
            --outputbasename {params.bin_prefix} \
            --search_engine {params.search_engine} \
            --write_bin_evals \
            --write_bins \
            --write_unbinned \
            --score_threshold {params.score_threshold} \
            --duplicate_penalty {params.duplicate_penalty} \
            --megabin_penalty {params.megabin_penalty} \
            --threads {threads} --debug > {log} 2>&1

            rm -rf ${{FNA%.gz}}
            rm -rf ${{PEP%.gz}}

            exitcode=$?
            if [ $exitcode -eq 1 ]
            then
                grep -oEi 'no single copy genes found. Aborting' {log}
                grepcode=$?
                if [ $grepcode -eq 0 ]
                then
                    exit 0
                else
                    grep -oEi 'single copy gene prediction using {params.search_engine} failed. Aborting' {log}
                    grepcode=$?
                    if [ $grepcode -eq 0 ]
                    then
                        exit 0
                    else
                        exit $exitcode
                    fi
                fi
            fi

            python {params.wrapper_dir}/dastools_postprocess.py \
            {params.bin_prefix}

            for FILESTR in `ls {params.bin_prefix}*`
            do
                if [ -f $FILESTR ] && [ "${{FILESTR##*.}}" != "gz" ];
                then
                    pigz -f $FILESTR
                fi
            done

            touch {output}
            ''' 


    rule binning_dastools_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "mags/{binning_group}.{assembly_group}.{assembler}/dastools/binning_done"),
                    zip,
                    binning_group=ASSEMBLY_GROUPS["binning_group"],
                    assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                    assembler=ASSEMBLY_GROUPS["assembler"])

            #rules.predict_scaftigs_gene_prodigal_all.input

else:
    rule binning_dastools_all:
        input:


localrules:
    binning_graphbin2_all,
    binning_dastools_all