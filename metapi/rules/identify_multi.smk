# reference: https://github.com/RasmussenLab/phamb/blob/master/workflows/mag_annotation/Snakefile 

if config["params"]["identify"]["phamb"]["do"]:
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
            min_contig = config["params"]["binning"]["vamb"]["min_contig"],
            hmmsearch_evalue = config["params"]["identify"]["phamb"]["hmmsearch_evalue"]
        threads:
            config["params"]["identify"]["threads"]
        conda:
            config["envs"]["phamb"]
        shell:
            '''
            hmmsearch \
            --cpu {threads} \
            -E {params.hmmsearch_evalue} \
            -o {params.tmp} \
            --tblout {output.hmm} \
            {input.db} \
            {input.pep} \
            2> {log}
            '''


    rule identify_phamb_micomplete_merge:
        input:
            lambda wildcards: expand(
                os.path.join(
                    config["output"]["identify"],
                    "phamb/hmm/{assembly_group}.{{assembler}}.hmmsearch.out/{assembly_group}.{{assembler}}.hmmMiComplete105.tbl"),
                    assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            os.path.join(config["output"]["identify"], "phamb/hmm_merge/{binning_group}.{assembler}.hmmMiComplete105.tbl")
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
            min_contig = config["params"]["binning"]["vamb"]["min_contig"],
            hmmsearch_evalue = config["params"]["identify"]["phamb"]["hmmsearch_evalue"]
        threads:
            config["params"]["identify"]["threads"]
        conda:
            config["envs"]["phamb"]
        shell:
            '''
            hmmsearch \
            --cpu {threads} \
            -E {params.hmmsearch_evalue} \
            -o {params.tmp} \
            --tblout {output.hmm} \
            {input.db} \
            {input.pep} \
            2> {log}
            '''


    rule identify_phamb_vog_merge:
        input:
            lambda wildcards: expand(os.path.join(
                config["output"]["identify"],
                "phamb/hmm/{assembly_group}.{{assembler}}.hmmsearch.out/{assembly_group}.{{assembler}}.hmmVOG.tbl"),
                assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            os.path.join(config["output"]["identify"], "phamb/hmm_merge/{binning_group}.{assembler}.hmmVOG.tbl")
        log:
            os.path.join(config["output"]["identify"], "logs/phamb_vog_merge.{binning_group}.{assembler}.log")
        shell:
            '''
            cat {input} > {output} 2> {log}
            '''


    rule identify_phamb_all:
        input:
            expand(
                os.path.join(config["output"]["identify"],
                             "phamb/hmm_merge/{binning_group}.{assembler}.{suffix}.tbl"),
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                suffix=["hmmMiComplete105", "hmmVOG"])

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
    identify_phamb_filter_pep,
    identify_phamb_micomplete_merge,
    identify_phamb_vog_merge,
    identify_phamb_all,
    identify_multi_all,
    identify_all