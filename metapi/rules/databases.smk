rule databases_bacteriome_refine_taxonomy:
    input:
        genomes_info = os.path.join(
            config["output"]["check"],
            "report/checkm/checkm_table_genomes_info.{assembler}.all.tsv"),
        table_gtdb = os.path.join(
            config["output"]["taxonomic"],
            "report/gtdbtk/MAGs_hmq.rep.{assembler}.{dereper}.gtdbtk.gtdb.tsv")
    output:
        taxonomy = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxonomy.{assembler}.{dereper}/taxonomy.tsv"),
        fna = directory(os.path.join(config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}"))
    params:
        db_name = config["params"]["databases"]["bacteriome"]["name"],
        rep_level = config["params"]["databases"]["bacteriome"]["rep_level"],
        out_dir = os.path.join(config["output"]["databases"], "bacteriome"),
        base_dir = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}")
    run:
        import shutil

        shutil.rmtree(params.out_dir) 

        metapi.refine_taxonomy(
            input.genomes_info,
            input.table_gtdb,
            params.db_name,
            params.rep_level,
            params.base_dir,
            output.taxonomy)


rule databases_bacteriome_generate_taxdump:
    input:
        taxonomy = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxonomy.{assembler}.{dereper}/taxonomy.tsv")
    output:
        taxdump = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
            taxdump=["names.dmp", "nodes.dmp", "taxid.map", "prelim_map.txt"]),
        taxmap = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/data/{taxmap}.map"),
            taxmap=["taxid", "name"])
    log:
        os.path.join(config["output"]["databases"],
            "logs/taxonkit/taxonkit.{assembler}.{dereper}.log")
    benchmark:
        os.path.join(config["output"]["databases"],
            "benchmark/taxonkit/taxonkit.{assembler}.{dereper}.benchmark.txt")
    params:
        out_dir = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}")
    conda:
        config["envs"]["taxonkit"]
    shell:
        '''
        rm -rf {params.out_dir}

        cat {input.taxonomy} | \
        csvtk cut -t -f genome_id,kingdom,phylum,class,order,family,genus,species,strain | \
        csvtk del-header | \
        taxonkit create-taxdump - \
        --out-dir {params.out_dir} \
        --force -A 1 \
        -R superkingdom,phylum,class,order,family,genus,species,strain 2> {log}

        awk '{{print "TAXID\t" $0}}' {params.out_dir}/taxid.map > {params.out_dir}/prelim_map.txt

        mkdir -p {params.out_dir}/data

        cp {params.out_dir}/taxid.map {params.out_dir}/data/taxid.map

        cat {params.out_dir}/taxid.map | \
        taxonkit lineage \
        --data-dir {params.out_dir} \
        -i 2 -n -L | \
        cut -f 1,3 > {params.out_dir}/data/name.map 2>> {log}
        '''


rule databases_bacteriome_extract_taxonomy:
    input:
        taxdump = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
            taxdump=["names.dmp", "nodes.dmp", "taxid.map", "prelim_map.txt"])
    output:
        os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab")
    log:
        os.path.join(config["output"]["databases"],
            "logs/krona/krona_extract_taxonomy.{assembler}.{dereper}.log")
    benchmark:
        os.path.join(config["output"]["databases"],
            "benchmark/krona/krona_extract_taxonomy.{assembler}.{dereper}.benchmark.txt")
    params:
        taxdumpdir = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}")
    conda:
        config["envs"]["kraken2"]
    shell:
        '''
        kronadir=$(dirname $(realpath $(which ktUpdateTaxonomy.sh)))
        extaxpl=$kronadir/scripts/extractTaxonomy.pl

        perl $extaxpl {params.taxdumpdir} >{log} 2>&1
        '''


rule databases_bacteriome_kmcp_compute:
    input:
        os.path.join(config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}")
    output:
        directory(os.path.join(config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kmcp/compute"))
    log:
        os.path.join(config["output"]["databases"],
            "logs/kmcp/kmcp_compute.{assembler}.{dereper}.log")
    benchmark:
        os.path.join(config["output"]["databases"],
            "benchmark/kmcp/kmcp_compute.{assembler}.{dereper}.benchmark.txt")
    params:
        kmer = config["params"]["databases"]["bacteriome"]["kmcp"]["compute"]["kmer"],
        split_number = config["params"]["databases"]["bacteriome"]["kmcp"]["compute"]["split_number"]
    conda:
        config["envs"]["kmcp"]
    threads:
        config["params"]["databases"]["threads"]
    shell:
        '''
        rm -rf {output}

        kmcp compute \
        --threads {threads} \
        --in-dir {input}/ \
        --out-dir {output} \
        --kmer {params.kmer} \
        --split-number {params.split_number} \
        --force >{log} 2>&1
        '''


rule databases_bacteriome_kmcp_index:
    input:
        taxmap = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/data/{taxmap}.map"),
            taxmap=["taxid", "name"]),
        fna = os.path.join(config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kmcp/compute")
    output:
        directory(os.path.join(config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kmcp/index"))
    log:
        os.path.join(config["output"]["databases"],
            "logs/kmcp/kmcp_index.{assembler}.{dereper}.log")
    benchmark:
        os.path.join(config["output"]["databases"],
            "benchmark/kmcp/kmcp_index.{assembler}.{dereper}.benchmark.txt")
    params:
        false_positive_rate = config["params"]["databases"]["bacteriome"]["kmcp"]["index"]["false_positive_rate"]
    conda:
        config["envs"]["kmcp"]
    threads:
        config["params"]["databases"]["threads"]
    shell:
        '''
        rm -rf {output}

        kmcp index \
        --threads {threads} \
        --in-dir {input.fna} \
        --out-dir {output} \
        --false-positive-rate {params.false_positive_rate} \
        --log {log}

        cp {input.taxmap} {output}/
        '''


rule databases_all:
    input:
        expand([
            os.path.join(
                config["output"]["databases"],
                "bacteriome/taxonomy.{assembler}.{dereper}/taxonomy.tsv"),
            os.path.join(
                config["output"]["databases"],
                "bacteriome/fna.{assembler}.{dereper}"),
            os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}/{taxdump}"),
            os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}/data/{taxmap}.map"),
            os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab"),
            os.path.join(
                config["output"]["databases"],
                "bacteriome/databases.{assembler}.{dereper}/kmcp/compute"),
            os.path.join(
                config["output"]["databases"],
                "bacteriome/databases.{assembler}.{dereper}/kmcp/index")],
            taxdump=["names.dmp", "nodes.dmp", "taxid.map", "prelim_map.txt"],
            taxmap=["taxid", "name"],
            assembler=ASSEMBLERS,
            dereper=DEREPERS)