# reference: https://github.com/RasmussenLab/phamb/blob/master/workflows/mag_annotation/Snakefile 
rule identify_phamb_filter_pep:
    input:
        pep = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binning_group}.{assembly_group}.{assembler}.prodigal.faa.gz"),
        gff = os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binning_group}.{assembly_group}.{assembler}.prodigal.gff.gz")
    output:
        pep = expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene_ge{min_contig}/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal.renamed.ge{min_contig}.faa.gz"),
            min_contig=config["params"]["binning"]["vamb"]["min_contig"])
    params:
        min_contig = config["params"]["binning"]["vamb"]["min_contig"],
        binning_group = "{binning_group}",
        assembly_group = "{assembly_group}"
    run:
        assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))
        assembly_index = int(assembly_groups.index(params.assembly_group)) + 1
        assembly_group = f'''S{assembly_index}'''

        pep_id_list = metapi.parse_gff(str(input.gff), int(params.min_contig))
        metapi.extract_faa(str(input.pep), pep_id_list, str(output.pep), assembly_group)


rule identify_phamb_hmmsearch_micomplete:
    input:
        pep = expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene_ge{min_contig}/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal.renamed.ge{min_contig}.faa.gz"),
            min_contig=config["params"]["binning"]["vamb"]["min_contig"]),
        db = config["params"]["identify"]["phamb"]["micompletedb"]
    output:
        hmm = os.path.join(
            config["output"]["identify"],
            "hmm/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.hmmMiComplete105.tbl.gz")
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_phamb_hmmsearch_micomplete/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_phamb_hmmsearch_micomplete/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        hmmsearch_evalue = config["params"]["identify"]["phamb"]["hmmsearch_evalue"],
        tmp = os.path.join(
            config["output"]["identify"],
            "hmm/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.hmmMiComplete105.tmp")
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["phamb"]
    shell:
        '''
        set +e

        HMMGZ={output.hmm}
        HMM=${{HMMGZ%.gz}}
        PEPGZ={input.pep}
        PEP=${{PEPGZ%.gz}}

        if [[ `zcat {input.pep} | wc -l` -eq 0 ]];
        then
            touch $HMM >> {log} 2>&1
            pigz -f $HMM >> {log} 2>&1
            exit 0
        else
            pigz -dkf $PEPGZ

            hmmsearch \
            --cpu {threads} \
            -E {params.hmmsearch_evalue} \
            -o {params.tmp} \
            --tblout $HMM \
            {input.db} \
            $PEP \
            >> {log} 2>&1

            rm -rf {params.tmp}
            rm -rf $PEP

            exitcode=$?
            if [ $exitcode -eq 0 ];
            then
                grep -oEi "ok" $HMM
                grepcode=$?
                if [ $grepcode -eq 1 ];
                then
                    echo "hmmsearch failed" >> {log} 2>&1
                    exit 1
                else
                    echo "hmmsearch done" >> {log} 2>&1
                    pigz -f $HMM
                    exit 0
                fi
            else
                echo "hmmsearch failed" >> {log} 2>&1
                exit 1
            fi
        fi
        '''


rule identify_phamb_hmmsearch_micomplete_merge:
    input:
        lambda wildcards: expand(
            os.path.join(
                config["output"]["identify"],
                "hmm/{{binning_group}}.{assembly_group}.{{assembler}}/{{binning_group}}.{assembly_group}.{{assembler}}.hmmMiComplete105.tbl.gz"),
                assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
    output:
        os.path.join(config["output"]["identify"], "annotations/{binning_group}.{assembler}/all.hmmMiComplete105.tbl.gz")
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_phamb_hmmsearch_micomplete_merge/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_phamb_hmmsearch_micomplete_merge/{binning_group}.{assembler}.txt")
    shell:
        '''
        cat {input} >{output} 2> {log}
        '''


rule identify_phamb_hmmsearch_vog:
    input:
        pep = expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene_ge{min_contig}/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal/{{binning_group}}.{{assembly_group}}.{{assembler}}.prodigal.renamed.ge{min_contig}.faa.gz"),
            min_contig=config["params"]["binning"]["vamb"]["min_contig"]),
        db = config["params"]["identify"]["phamb"]["vogdb"]
    output:
        hmm = os.path.join(
            config["output"]["identify"],
            "hmm/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.hmmVOG.tbl.gz")
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_phamb_hmmsearch_vog/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_phamb_hmmsearch_vog/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        hmmsearch_evalue = config["params"]["identify"]["phamb"]["hmmsearch_evalue"],
        tmp = os.path.join(
            config["output"]["identify"],
            "hmm/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.hmmVOG.tmp")
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["phamb"]
    shell:
        '''
        set +e

        HMMGZ={output.hmm}
        HMM=${{HMMGZ%.gz}}
        PEPGZ={input.pep}
        PEP=${{PEPGZ%.gz}}

        if [[ `zcat {input.pep} | wc -l` -eq 0 ]];
        then
            touch $HMM >> {log} 2>&1
            pigz -f $HMM >> {log} 2>&1
            exit 0
        else
            pigz -dkf $PEPGZ

            hmmsearch \
            --cpu {threads} \
            -E {params.hmmsearch_evalue} \
            -o {params.tmp} \
            --tblout $HMM \
            {input.db} \
            $PEP \
            >> {log} 2>&1

            rm -rf {params.tmp}
            rm -rf $PEP

            exitcode=$?
            if [ $exitcode -eq 0 ];
            then
                grep -oEi "ok" $HMM
                grepcode=$?
                if [ $grepcode -eq 1 ];
                then
                    echo "hmmsearch failed" >> {log} 2>&1
                    exit 1
                else
                    echo "hmmsearch done" >> {log} 2>&1
                    pigz -f $HMM
                    exit 0
                fi
            else
                echo "hmmsearch failed" >> {log} 2>&1
                exit 1
            fi
        fi
        '''


rule identify_phamb_hmmsearch_vog_merge:
    input:
        lambda wildcards: expand(os.path.join(
            config["output"]["identify"],
            "hmm/{{binning_group}}.{assembly_group}.{{assembler}}/{{binning_group}}.{assembly_group}.{{assembler}}.hmmVOG.tbl.gz"),
            assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
    output:
        os.path.join(config["output"]["identify"], "annotations/{binning_group}.{assembler}/all.hmmVOG.tbl.gz")
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_phamb_hmmsearch_vog_merge/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_phamb_hmmsearch_vog_merge/{binning_group}.{assembler}.txt")
    shell:
        '''
        cat {input} >{output} 2> {log}
        '''


rule identify_phamb_randomforest:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}/binning_done"),
        micomplete = os.path.join(
            config["output"]["identify"],
            "annotations/{binning_group}.{assembler}/all.hmmMiComplete105.tbl.gz"),
        vog = os.path.join(
            config["output"]["identify"],
            "annotations/{binning_group}.{assembler}/all.hmmVOG.tbl.gz"),
        dvf = os.path.join(
            config["output"]["identify"],
            "annotations/{binning_group}.{assembler}/all.DVF.predictions.txt.gz")
    output:
        # vambbins_aggregated_annotation.txt
        # vambbins_RF_predictions.txt
        # vamb_bins
        os.path.join(config["output"]["identify"], "vmags_phamb/{binning_group}.{assembler}/phamb_randomforest_done")
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_phamb_randomforest/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_phamb_randomforest/{binning_group}.{assembler}.txt")
    params:
        randomforest_script = config["params"]["identify"]["phamb"]["randomforest_script"],
        min_binsize = config["params"]["identify"]["phamb"]["min_binsize"],
        mags_dir = os.path.join(config["output"]["binning"], "mags_vamb/{binning_group}.{assembler}"),
        annotations_dir = os.path.join(config["output"]["identify"], "annotations/{binning_group}.{assembler}"),
        output_dir = os.path.join(config["output"]["identify"], "vmags_phamb/{binning_group}.{assembler}")
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["phamb"]
    shell:
        '''
        if [ -e {params.mags_dir}/clusters.tsv ];
        then
            python {params.randomforest_script} \
            -m {params.min_binsize} \
            -s C \
            {input.scaftigs} \
            {params.mags_dir}/clusters.tsv \
            {params.annotations_dir} \
            {params.output_dir} \
            >{log} 2>&1

            touch {output}
        else
            touch {output}
        fi
        '''


rule identify_phamb_postprocess:
    input:
        metadata = os.path.join(config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz"),
        phamb_rf_done = os.path.join(config["output"]["identify"],
            "vmags_phamb/{binning_group}.{assembler}/phamb_randomforest_done")
    output:
        viral = os.path.join(
            config["output"]["identify"],
            "vmags/{binning_group}.{assembly_group}.{assembler}/phamb/{binning_group}.{assembly_group}.{assembler}.phamb.combined.fa.gz")
    params:
        binning_group = "{binning_group}",
        assembly_group = "{assembly_group}",
        assembler = "{assembler}"
    run:
        import os
        import gzip
        import subprocess
        from glob import glob
        from Bio import SeqIO

        binning_assembly_metadata = pd.read_csv(input.metadata, sep="\t").set_index("binning_assembly_group")
        assembly_index = binning_assembly_metadata.loc[f'''{params.binning_group}.{params.assembly_group}''', "vamb_id"]

        vamb_bins_dir = os.path.join(os.path.dirname(input.phamb_rf_done), "vamb_bins")

        os.makedirs(os.path.dirname(output.viral), exist_ok=True)

        if os.path.exists(vamb_bins_dir):
            with gzip.open(output.viral, 'wt') as oh:
                vamb_bins_list = glob(f'''{vamb_bins_dir}/vamb_bins.*.fna''')
                if len(vamb_bins_list) > 0:
                    for fna in vamb_bins_list:
                        for rc in SeqIO.parse(fna, "fasta"):
                            if rc.id.startswith(f'''{assembly_index}C'''):
                                SeqIO.write(rc, oh, "fasta")
                else:
                    subprocess.run(f'''touch {output.viral}''', shell=True)
        else:
            subprocess.run(f'''touch {output.viral}''', shell=True)


if config["params"]["identify"]["phamb"]["do"] and \
config["params"]["identify"]["deepvirfinder"]["do"] and \
config["params"]["binning"]["vamb"]["do"]:
    rule identify_phamb_all:
        input:
            expand([
                os.path.join(
                    config["output"]["identify"],
                    "annotations/{binning_group}.{assembler}/all.{annotations}.gz"),
                os.path.join(
                    config["output"]["identify"],
                    "vmags_phamb/{binning_group}.{assembler}/phamb_randomforest_done")],
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                annotations=["hmmMiComplete105.tbl", "hmmVOG.tbl"]),
            expand(
                os.path.join(
                    config["output"]["identify"],
                    "vmags/{binning_group}.{assembly_group}.{assembler}/phamb/{binning_group}.{assembly_group}.{assembler}.phamb.combined.fa.gz"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"])

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
