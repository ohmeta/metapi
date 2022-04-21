# reference: https://github.com/RasmussenLab/phamb/blob/master/workflows/mag_annotation/Snakefile 

if config["params"]["identify"]["phamb"]["do"]:
    rule identify_phamb_filter_pep:
        input:
            pep = os.path.join(config["output"]["predict"],
                               "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.faa"),
            gff = os.path.join(config["output"]["predict"],
                               "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.gff")
        output:   
            pep = expand(os.path.join(config["output"]["predict"],
                         "scaftigs_gene/{assembly_group}.{assembler}.prodigal.out/{assembly_group}.{assembler}.faa.ge{min_contig}"),
                         min_contig=config["params"]["binning"]["vamb"]["min_contig"])
        params:
            min_contig = config["params"]["binning"]["vamb"]["min_contig"]
        run:
            pep_id_list = metapi.parse(input.gff, int(params.min_contig))
            metapi.extract_faa(input.pep, pep_id_list, output.pep)
      

    rule identify_phamb_micomplete:
        input:
            pep = expand(os.path.join(config["output"]["predict"],
                         "scaftigs_gene/{{assembly_group}}.{{assembler}}.prodigal.out/{{assembly_group}}.{{assembler}}.faa.ge{min_contig}"),
                         min_contig=config["params"]["binning"]["vamb"]["min_contig"]),
            db = config["params"]["classify"]["phamb"]["micompletedb"]
        output:
            hmm = os.path.join(config["output"]["classify"],
                               "hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.hmmMiComplete105.tbl")
        log:
            os.path.join(config["output"]["classify"],
                         "logs/phamb_hmmsearch_micomplete/{assembly_group}.{assembler}.phamb_hmmsearch_micomplete.log")
        benchmark:
            os.path.join(config["output"]["classify"],
                         "benchmark/phamb_hmmsearch_micomplete/{assembly_group}.{assembler}.phamb_hmmsearch_micomplete.log")
        params:
            tmp = os.path.join(config["output"]["classify"],
                               "hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.micomplete.tmp"),
            min_contig = config["params"]["binning"]["vamb"]["min_contig"],
            hmmsearch_evalue = config["params"]["classify"]["phamb"]["hmmsearch_evalue"]
        threads:
            config["params"]["classify"]["threads"]
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


    rule identify_phamb_vog:
        input:
            pep = expand(os.path.join(config["output"]["predict"],
                         "scaftigs_gene/{{assembly_group}}.{{assembler}}.prodigal.out/{{assembly_group}}.{{assembler}}.faa.ge{min_contig}"),
                         min_contig=config["params"]["binning"]["vamb"]["min_contig"]),
            db = config["params"]["classify"]["phamb"]["vogdb"]
        output:
            hmm = os.path.join(config["output"]["classify"],
                               "hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.hmmVOG.tbl")
        log:
            os.path.join(config["output"]["classify"],
                         "logs/phamb_hmmsearch_vog/{assembly_group}.{assembler}.phamb_hmmsearch_vog.log")
        benchmark:
            os.path.join(config["output"]["classify"],
                         "benchmark/phamb_hmmsearch_vog/{assembly_group}.{assembler}.phamb_hmmsearch_vog.log")
        params:
            tmp = os.path.join(config["output"]["classify"],
                               "hmm/{assembly_group}.{assembler}.hmmsearch.out/{assembly_group}.{assembler}.micomplete.tmp"),
            min_contig = config["params"]["binning"]["vamb"]["min_contig"],
            hmmsearch_evalue = config["params"]["classify"]["phamb"]["hmmsearch_evalue"]
        threads:
            config["params"]["classify"]["threads"]
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