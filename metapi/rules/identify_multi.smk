# reference: https://github.com/RasmussenLab/phamb/blob/master/workflows/mag_annotation/Snakefile 

if config["params"]["identify"]["phamb"]["do"] and config["params"]["identify"]["deepvirfinder"]["do"] and config["params"]["binning"]["vamb"]["do"]:
    rule identify_phamb_filter_pep:
        input:
            pep = os.path.join(config["output"]["predict"],
                               "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.faa"),
            gff = os.path.join(config["output"]["predict"],
                               "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.gff")
        output:
            pep = expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{assembly_group}}.{{assembler}}.prodigal.out/{{assembly_group}}.{{assembler}}.renamed.ge{min_contig}.faa"),
                         min_contig=config["params"]["binning"]["vamb"]["min_contig"])
        params:
            min_contig = config["params"]["binning"]["vamb"]["min_contig"],
            assembly_group = "{assembly_group}"
        run:
            binning_group = metapi.get_binning_group_by_assembly_group(SAMPLES, params.assembly_group)[0]
            assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, binning_group))
            assembly_index = int(assembly_groups.index(params.assembly_group)) + 1
            assembly_group = f'''S{assembly_index}'''

            pep_id_list = metapi.parse_gff(str(input.gff), int(params.min_contig))
            metapi.extract_faa(str(input.pep), pep_id_list, str(output.pep), assembly_group) 


    rule identify_phamb_micomplete:
        input:
            pep = expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{assembly_group}}.{{assembler}}.prodigal.out/{{assembly_group}}.{{assembler}}.renamed.ge{min_contig}.faa"),
                         min_contig=config["params"]["binning"]["vamb"]["min_contig"]),
            db = config["params"]["identify"]["phamb"]["micompletedb"]
        output:
            hmm = os.path.join(config["output"]["identify"],
                               "phamb/hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.hmmMiComplete105.tbl")
        log:
            os.path.join(config["output"]["identify"],
                         "logs/phamb_hmmsearch_micomplete/{assembly_group}.{assembler}.phamb_hmmsearch_micomplete.log")
        benchmark:
            os.path.join(config["output"]["identify"],
                         "benchmark/phamb_hmmsearch_micomplete/{assembly_group}.{assembler}.phamb_hmmsearch_micomplete.benchmark.txt")
        params:
            tmp = os.path.join(config["output"]["identify"],
                               "phamb/hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.hmmMiComplete105.tmp"),
            hmmsearch_evalue = config["params"]["identify"]["phamb"]["hmmsearch_evalue"]
        threads:
            config["params"]["identify"]["threads"]
        conda:
            config["envs"]["phamb"]
        shell:
            '''
            set +e
            if [[ `cat {input.pep} | wc -l` -eq 0 ]];
            then
                touch {output.hmm} >> {log} 2>&1
                exit 0
            else
                hmmsearch \
                --cpu {threads} \
                -E {params.hmmsearch_evalue} \
                -o {params.tmp} \
                --tblout {output.hmm} \
                {input.db} \
                {input.pep} \
                >> {log} 2>&1

                exitcode=$?
                if [ $exitcode -eq 0 ];
                then
                    grep -oEi "ok" {output.hmm}
                    grepcode=$?
                    if [ $grepcode -eq 1 ];
                    then
                        rm -rf {output.hmm}
                        echo "hmmsearch failed" >> {log} 2>&1
                        exit 1
                    else
                        echo "hmmsearch done" >> {log} 2>&1
                        exit 0
                    fi
                else
                    rm -rf {output.hmm}
                    echo "hmmsearch failed" >> {log} 2>&1
                    exit 1
                fi
            fi
            '''


    rule identify_phamb_micomplete_merge:
        input:
            lambda wildcards: expand(
                os.path.join(
                    config["output"]["identify"],
                    "phamb/hmm/{assembly_group}.{{assembler}}.hmmsearch.out/{assembly_group}.{{assembler}}.hmmMiComplete105.tbl"),
                    assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            os.path.join(config["output"]["identify"], "phamb/annotations/{binning_group}.{assembler}/all.hmmMiComplete105.tbl")
        log:
            os.path.join(config["output"]["identify"], "logs/phamb_micomplete_merge.{binning_group}.{assembler}.log")
        shell:
            '''
            cat {input} > {output} 2> {log}
            '''
 

    rule identify_phamb_vog:
        input:
            pep = expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{assembly_group}}.{{assembler}}.prodigal.out/{{assembly_group}}.{{assembler}}.renamed.ge{min_contig}.faa"),
                         min_contig=config["params"]["binning"]["vamb"]["min_contig"]),
            db = config["params"]["identify"]["phamb"]["vogdb"]
        output:
            hmm = os.path.join(config["output"]["identify"],
                               "phamb/hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.hmmVOG.tbl")
        log:
            os.path.join(config["output"]["identify"],
                         "logs/phamb_hmmsearch_vog/{assembly_group}.{assembler}.phamb_hmmsearch_vog.log")
        benchmark:
            os.path.join(config["output"]["identify"],
                         "benchmark/phamb_hmmsearch_vog/{assembly_group}.{assembler}.phamb_hmmsearch_vog.benchmark.txt")
        params:
            tmp = os.path.join(config["output"]["identify"],
                               "phamb/hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.hmmVOG.tmp"),
            hmmsearch_evalue = config["params"]["identify"]["phamb"]["hmmsearch_evalue"]
        threads:
            config["params"]["identify"]["threads"]
        conda:
            config["envs"]["phamb"]
        shell:
            '''
            if [[ `cat {input.pep} | wc -l` -eq 0 ]];
            then
                touch {output.hmm} >> {log} 2>&1
                exit 0
            else
                hmmsearch \
                --cpu {threads} \
                -E {params.hmmsearch_evalue} \
                -o {params.tmp} \
                --tblout {output.hmm} \
                {input.db} \
                {input.pep} \
                >> {log} 2>&1

                exitcode=$?
                if [ $exitcode -eq 0 ];
                then
                    grep -oEi "ok" {output.hmm}
                    grepcode=$?
                    if [ $grepcode -eq 1 ];
                    then
                        rm -rf {output.hmm}
                        echo "hmmsearch failed" >> {log} 2>&1
                        exit 1
                    else
                        echo "hmmsearch done" >> {log} 2>&1
                        exit 0
                    fi
                else
                    rm -rf {output.hmm}
                    echo "hmmsearch failed" >> {log} 2>&1
                    exit 1
                fi
            fi
            '''


    rule identify_phamb_vog_merge:
        input:
            lambda wildcards: expand(os.path.join(
                config["output"]["identify"],
                "phamb/hmm/{assembly_group}.{{assembler}}.hmmsearch.out/{assembly_group}.{{assembler}}.hmmVOG.tbl"),
                assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            os.path.join(config["output"]["identify"], "phamb/annotations/{binning_group}.{assembler}/all.hmmVOG.tbl")
        log:
            os.path.join(config["output"]["identify"], "logs/phamb_vog_merge.{binning_group}.{assembler}.log")
        shell:
            '''
            cat {input} > {output} 2> {log}
            '''


    rule identify_phamb_deepvirfinder_merge:
        input:
            lambda wildcards: expand(
                os.path.join(
                    config["output"]["identify"],
                    "deepvirfinder/{assembly_group}.{{assembler}}.dvf.out/{assembly_group}.{{assembler}}.scaftigs.fa.gz_gt{min_length}bp_dvfpred.txt"),
                min_length=config["params"]["identify"]["deepvirfinder"]["min_length"],
                assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            dvf = os.path.join(config["output"]["identify"], "phamb/annotations/{binning_group}.{assembler}/all.DVF.predictions.txt")
        params:
            binning_group = "{binning_group}"
        run:
            assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))
            with open(output.dvf, "w") as oh:
                assembly_group = os.path.basename(input[0]).split(".")[0]
                assembly_index = int(assembly_groups.index(assembly_group)) + 1
                assembly_group = f'''S{assembly_index}'''
                with open(input[0], "r") as ih:
                    oh.write(ih.readline())
                    for line in ih:
                        oh.write(f'''{assembly_group}C{line}''')
                for dvfpred in input[1:]:
                    assembly_group = os.path.basename(dvfpred).split(".")[0]
                    assembly_index = int(assembly_groups.index(assembly_group)) + 1
                    assembly_group = f'''S{assembly_index}'''
                    with open(dvfpred, "r") as ih:
                        ih.readline()
                        for line in ih:
                            oh.write(f'''{assembly_group}C{line}''')


    rule identify_phamb_randomforest:
        input:
            scaftigs = os.path.join(
                config["output"]["multisplit_binning"],
                "scaftigs/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.fa.gz"),
            binning_done = os.path.join(
                config["output"]["multisplit_binning"],
                "bins/{binning_group}.{assembler}.vamb.out/binning_done"),
            micomplete = os.path.join(
                config["output"]["identify"],
                "phamb/annotations/{binning_group}.{assembler}/all.hmmMiComplete105.tbl"),
            vog = os.path.join(
                config["output"]["identify"],
                "phamb/annotations/{binning_group}.{assembler}/all.hmmVOG.tbl"),
            dvf = os.path.join(
                config["output"]["identify"],
                "phamb/annotations/{binning_group}.{assembler}/all.DVF.predictions.txt")
        output:
            #annotation = os.path.join(config["output"]["identify"], "phamb/randomforest/{binning_group}.{assembler}/vambbins_aggregated_annotation.txt"),
            #prediction = os.path.join(config["output"]["identify"], "phamb/randomforest/{binning_group}.{assembler}/vambbins_RF_predictions.txt"),
            #vamb_bins = directory(os.path.join(config["output"]["identify"], "phamb/randomforest/{binning_group}.{assembler}/vamb_bins")),
            os.path.join(config["output"]["identify"], "phamb/randomforest/{binning_group}.{assembler}/phamb_randomforest_done")
        log:
            os.path.join(config["output"]["identify"],
                         "logs/phamb_randomforest/{binning_group}.{assembler}.phamb_randomforest.log")
        benchmark:
            os.path.join(config["output"]["identify"],
                         "benchmark/phamb_randomforest/{binning_group}.{assembler}.phamb_randomforest.benchmark.txt")
        conda:
            config["envs"]["phamb"]
        params:
            randomforest_script = config["params"]["identify"]["phamb"]["randomforest_script"],
            min_binsize = config["params"]["identify"]["phamb"]["min_binsize"],
            bins_dir = os.path.join(config["output"]["multisplit_binning"], "bins/{binning_group}.{assembler}.vamb.out"),
            annotations_dir = os.path.join(config["output"]["identify"], "phamb/annotations/{binning_group}.{assembler}"),
            output_dir = os.path.join(config["output"]["identify"], "phamb/randomforest/{binning_group}.{assembler}")
        shell:
            '''
            if [ -e {params.bins_dir}/clusters.tsv ];
            then
                python {params.randomforest_script} \
                -m {params.min_binsize} \
                -s C \
                {input.scaftigs} \
                {params.bins_dir}/clusters.tsv \
                {params.annotations_dir} \
                {params.output_dir} \
                > {log} 2>&1

                touch {output}
            else
                touch {output}
            fi
            '''


    rule identify_phamb_all:
        input:
            expand([
                os.path.join(config["output"]["identify"],
                             "phamb/annotations/{binning_group}.{assembler}/all.{annotations}"),
                os.path.join(config["output"]["identify"],
                             "phamb/randomforest/{binning_group}.{assembler}/phamb_randomforest_done")],
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                annotations=["hmmMiComplete105.tbl", "hmmVOG.tbl", "DVF.predictions.txt"])
            
else:
    rule identify_phamb_all:
        input:


rule identify_multi_all:
    input:
        rules.identify_phamb_all.input


rule identify_all:
    input:
        rules.identify_single_all.input,
        rules.identify_multi_all.input


localrules:
    identify_phamb_all,
    identify_multi_all,
    identify_all