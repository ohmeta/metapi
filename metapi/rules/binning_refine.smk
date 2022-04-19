if config["params"]["binning"]["graphbin2"]["do"]:
    rule binning_graphbin2_prepare_assembly:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz"),
            gfa = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.gfa.gz"),
        output:
             scaftigs = temp(os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/graphbin2/scaftigs.fa")),
             gfa = temp(os.path.join(
                 config["output"]["binning"],
                 "bins/{assembly_group}.{assembler}.out/graphbin2/scaftigs.gfa"))
        shell:
            '''
            pigz -dc {input.scaftigs} > {output.scaftigs}
            pigz -dc {input.gfa} > {output.gfa}
            '''

           
    rule binning_graphbin2_prepare_binned:
        input:
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner_graphbin}")
        output:
            binned = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/graphbin2/{assembly_group}.{assembler}.{binner_graphbin}.graphbin2.csv")
        params:
            assembler = "{assembler}"
        run:
            metapi.get_binning_info(input.bins_dir,
                                    output.binned,
                                    params.assembler)


    localrules:
        binning_graphbin2_prepare_binned


    rule binning_graphbin2:
        input:
            scaftigs = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/graphbin2/scaftigs.fa"),
            gfa = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/graphbin2/scaftigs.gfa"),
            binned = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/graphbin2/{assembly_group}.{assembler}.{binner_graphbin}.graphbin2.csv")
        output:
            os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner_graphbin}_graphbin2/binning_done")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{assembly_group}.{assembler}.{binner_graphbin}.graphbin2.refine.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/{binner_graphbin}/{assembly_group}.{assembler}.{binner_graphbin}.benchmark.txt")
        params:
            assembler = "{assembler}",
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner_graphbin}_graphbin2/"),
            prefix = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner_graphbin}_graphbin2/{assembly_group}.{assembler}.{binner_graphbin}_graphbin2.bin"),
            paths = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.paths.gz"),
            depth = config["params"]["binning"]["graphbin2"]["depth"],
            threshold = config["params"]["binning"]["graphbin2"]["threshold"]
        threads:
            config["params"]["binning"]["threads"]
        run:
            import pandas as pd
            import os

            shell('''mkdir -p {params.bins_dir}''')

            df = pd.read_csv(input.binned, names=["scaftigs_id", "bin_id"])

            if not df.empty:
                if params.assembler == "metaspades" or params.assembler == "spades":
                    shell(
                        '''
                        pigz -p {threads} -dc {params.paths} > {params.bins_dir}/scaftigs.paths

                        graphbin2 \
                        --assembler spades \
                        --contigs {input.scaftigs} \
                        --graph {input.gfa} \
                        --paths {params.bins_dir}/scaftigs.paths \
                        --binned {input.binned} \
                        --nthreads {threads} \
                        --depth {params.depth} \
                        --threshold {params.threshold} \
                        --output {params.bins_dir} \
                        > {log} 2>&1

                        rm -rf {params.bins_dir}/scaftigs.paths
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
                        --output {params.bins_dir} \
                        > {log} 2>&1
                        ''')

                metapi.generate_bins(f"{params.bins_dir}/graphbin2_output.csv",
                                     input.scaftigs, params.prefix)
                shell('''touch {output}''')


    rule binning_graphbin2_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner_graphbin}_graphbin2/binning_done"),
                binner_graphbin=BINNERS_GRAPHBIN,
                assembler=ASSEMBLERS,
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

            #rules.alignment_all.input,
            #rules.assembly_all.input

else:
    rule binning_graphbin2_all:
        input:


def get_binning_done(wildcards, binners):
    if len(binners) == 1:
        if binners[0] != "vamb":
            binning_done = expand(os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner}/binning_done"),
                assembly_group=wildcards.assembly_group,
                assembler=wildcards.assembler,
                binner=binners[0])
            return binning_done
        else:
            binning_done = expand(os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner}/binning_done"),
                assembly_group=metapi.get_multibinning_group_by_assembly_group(SAMPLES, wildcards.assembly_group),
                assembler=wildcards.assembler,
                binner=binners[0])
            return binning_done
    else:
        binning_done_list = []
        for binner in binners:
            if binner != "vamb":
                binning_done = expand(os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}.{assembler}.out/{binner}/binning_done"),
                    assembly_group=wildcards.assembly_group,
                    assembler=wildcards.assembler,
                    binner=binner)
                binning_done_list.append(binning_done[0])
            else:
                binning_done = expand(os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}.{assembler}.out/{binner}/binning_done"),
                    assembly_group=metapi.get_multibinning_group_by_assembly_group(SAMPLES, wildcards.assembly_group),
                    assembler=wildcards.assembler,
                    binner=binner)
                binning_done_list.append(binning_done[0])
        return binning_done_list
 

if config["params"]["binning"]["dastools"]["do"]:
    rule binning_dastools_preprocess:
        input:
            lambda wildcards: get_binning_done(wildcards, BINNERS_DASTOOLS)
        output:
            contigs2bin = expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins_id/{{assembly_group}}.{{assembler}}.out/{binner_dastools}_Contigs2Bin.tsv"),
                    binner_dastools=BINNERS_DASTOOLS)
        run:
            import glob
            import os
            from Bio import SeqIO

            i = -1
            for binning_done in input:
                i += 1
                bins_dir = os.path.dirname(binning_done)
                shell(f'''rm -rf {output.contigs2bin[i]}''')
                bins_list = glob.glob(bins_dir + "/*.bin.*.fa")
                if len(bins_list) == 0:
                    shell(f'''touch {output.contigs2bin[i]}''')
                else:
                    with open(output[i], 'w') as oh:
                        for bin_fa in sorted(bins_list):
                            bin_id_list = os.path.basename(bin_fa).split(".")
                            bin_id = bin_id_list[2] + "." + str(bin_id_list[4])
                            for contig in SeqIO.parse(bin_fa, "fasta"):
                                oh.write(f'''{contig.id}\t{bin_id}\n''')


    rule binning_dastools:
        input:
            contigs2bin = expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins_id/{{assembly_group}}.{{assembler}}.out/{binner_dastools}_Contigs2Bin.tsv"),
                    binner_dastools=BINNERS_DASTOOLS),
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz"),
            pep = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.faa")
        output:
            os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/dastools/binning_done")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{assembly_group}.{assembler}.dastools.binning.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/dastools/{assembly_group}.{assembler}.dastools.benchmark.txt")
        priority:
            30
        conda:
            config["envs"]["dastools"]
        params:
            binner = ",".join(BINNERS_DASTOOLS),
            wrapper_dir = WRAPPER_DIR,
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/dastools"),
            search_engine = config["params"]["binning"]["dastools"]["search_engine"],
            score_threshold = config["params"]["binning"]["dastools"]["score_threshold"],
            duplicate_penalty = config["params"]["binning"]["dastools"]["duplicate_penalty"],
            megabin_penalty = config["params"]["binning"]["dastools"]["megabin_penalty"],
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/dastools/{assembly_group}.{assembler}.dastools.bin")
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            set +e
            rm -rf {params.bins_dir}
            mkdir -p {params.bins_dir}
            
            pigz -p {threads} -d -c {input.scaftigs} > {params.bins_dir}/scaftigs.fasta
            contigs2bin=$(python -c "import sys; print(','.join(sys.argv[1:]))" {input.contigs2bin})

            DAS_Tool \
            --bins $contigs2bin \
            --labels {params.binner} \
            --contigs {params.bins_dir}/scaftigs.fasta \
            --proteins {input.pep} \
            --outputbasename {params.bin_prefix} \
            --search_engine {params.search_engine} \
            --write_bin_evals \
            --write_bins \
            --write_unbinned \
            --score_threshold {params.score_threshold} \
            --duplicate_penalty {params.duplicate_penalty} \
            --megabin_penalty {params.megabin_penalty} \
            --threads {threads} --debug > {log} 2>&1

            rm -rf {params.bins_dir}/scaftigs.fasta

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

            touch {output}
            ''' 


    rule binning_dastools_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}.{assembler}.out/dastools/binning_done"),
                assembler=ASSEMBLERS,
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

            #rules.predict_scaftigs_gene_prodigal_all.input

else:
    rule binning_dastools_all:
        input:


localrules:
    binning_graphbin2_all,
    binning_dastools_all