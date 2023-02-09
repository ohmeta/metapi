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
        fnadone = os.path.join(config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/done")
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

        shell("touch {output.fnadone}")


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
        taxtab = os.path.join(
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


if config["params"]["databases"]["bacteriome"]["kmcp"]["do"]:
    rule databases_bacteriome_kmcp_compute:
        input:
            fnadone = os.path.join(config["output"]["databases"],
                "bacteriome/fna.{assembler}.{dereper}/done")
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

            fnadir=$(dirname {input.fnadone})

            kmcp compute \
            --threads {threads} \
            --in-dir $fnadir/ \
            --out-dir {output} \
            --kmer {params.kmer} \
            --split-number {params.split_number} \
            --force >{log} 2>&1
            '''


    rule databases_bacteriome_kmcp_index:
        input:
            taxtab = os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab"),
            taxmap = expand(os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{{assembler}}.{{dereper}}/data/{taxmap}.map"),
                taxmap=["taxid", "name"]),
            compute = os.path.join(config["output"]["databases"],
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
            --in-dir {input.compute} \
            --out-dir {output} \
            --false-positive-rate {params.false_positive_rate} \
            --log {log}

            cp {input.taxmap} {output}/
            '''


    rule databases_bacteriome_kmcp_all:
        input:
            expand([
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/kmcp/compute"),
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/kmcp/index")],
            assembler=ASSEMBLERS,
            dereper=DEREPERS)

else:
    rule databases_bacteriome_kmcp_all:
        input:


if config["params"]["databases"]["bacteriome"]["kraken2"]["do"]:
    rule databases_bacteriome_kraken2_build:
        input:
            fnadone = os.path.join(config["output"]["databases"],
                "bacteriome/fna.{assembler}.{dereper}/done"),
            taxtab = os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab"),
            taxdump = expand(os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
                taxdump=["names.dmp", "nodes.dmp", "taxid.map", "prelim_map.txt"])
        output:
            expand(os.path.join(config["output"]["databases"],
                "bacteriome/databases.{{assembler}}.{{dereper}}/kraken2/index/{k2}.k2d"),
                k2=["hash", "opts", "taxo"])
        log:
            os.path.join(config["output"]["databases"],
                "logs/kraken2/kraken2_build.{assembler}.{dereper}.log")
        benchmark:
            os.path.join(config["output"]["databases"],
                "benchmark/kraken2/kraken2_build.{assembler}.{dereper}.benchmark.txt")
        params:
            tax = os.path.join(config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}"),
            db = os.path.join(config["output"]["databases"],
                "bacteriome/databases.{assembler}.{dereper}/kraken2/index")
        conda:
            config["envs"]["kraken2"]
        threads:
            config["params"]["databases"]["threads"]
        shell:
            '''
            dbdir=$(realpath {params.db})
            rm -rf $dbdir
            mkdir -p $dbdir

            fnadir=$(realpath $(dirname {input.fnadone}))
            dmpdir=$(realpath {params.tax})

            cp {input.taxdump[3]} $fnadir/

            pushd $dbdir
            ln -s $fnadir library
            ln -s $dmpdir taxonomy
            popd

            kraken2-build \
            --build \
            --threads {threads} \
            --db $dbdir \
            >{log} 2>&1
            '''


    rule databases_bacteriome_kraken2_bracken_build:
        input:
            expand(os.path.join(config["output"]["databases"],
                "bacteriome/databases.{{assembler}}.{{dereper}}/kraken2/index/{k2}.k2d"),
                k2=["hash", "opts", "taxo"])
        output:
            expand(os.path.join(config["output"]["databases"],
                "bacteriome/databases.{{assembler}}.{{dereper}}/kraken2/index/database{kmers}mers.{suffix}"),
                kmers=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["kmers"],
                suffix=["kmer_distrib", "kraken"])
        log:
            os.path.join(config["output"]["databases"],
                "logs/kraken2/bracken_build.{assembler}.{dereper}.log")
        benchmark:
            os.path.join(config["output"]["databases"],
                "benchmark/kraken2/bracken_build.{assembler}.{dereper}.benchmark.txt")
        params:
            db = os.path.join(config["output"]["databases"],
                "bacteriome/databases.{assembler}.{dereper}/kraken2/index"),
            ksize=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["ksize"],
            kmers=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["kmers"]
        conda:
            config["envs"]["kraken2"]
        threads:
            config["params"]["databases"]["threads"]
        shell:
            '''
            dbdir=$(realpath {params.db})

            rm -rf {log}

            for kmer in {params.kmers};
            do
                bracken-build -d $dbdir -t {threads} -k {params.ksize} -l $kmer >>{log} 2>&1
            done
            '''


    rule databases_bacteriome_kraken2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/kraken2/index/{k2}.k2d"),
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/kraken2/index/database{kmers}mers.{suffix}")],
                    k2=["hash", "opts", "taxo"],
                    kmers=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["kmers"],
                    suffix=["kmer_distrib", "kraken"],
                    assembler=ASSEMBLERS,
                    dereper=DEREPERS)

else:
    rule databases_bacteriome_kraken2_all:
        input:


if config["params"]["databases"]["bacteriome"]["krakenuniq"]["do"]:
    rule databases_bacteriome_krakenuniq_build:
        input:
            fnadone = os.path.join(config["output"]["databases"],
                "bacteriome/fna.{assembler}.{dereper}/done"),
            taxtab = os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab"),
            taxdump = expand(os.path.join(
                config["output"]["databases"],
                "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
                taxdump=["names.dmp", "nodes.dmp", "taxid.map", "prelim_map.txt"])
        output:
            expand(os.path.join(config["output"]["databases"],
                "bacteriome/databases.{{assembler}}.{{dereper}}/krakenuniq/index/{ku}"),
                ku=["database.kdb", "database.idx", "taxDB"])
        log:
            os.path.join(config["output"]["databases"],
                "logs/krakenuniq/krakenuniq_build.{assembler}.{dereper}.log")
        benchmark:
            os.path.join(config["output"]["databases"],
                "benchmark/krakenuniq/krakenuniq_build.{assembler}.{dereper}.benchmark.txt")
        params:
            tax = os.path.join(config["output"]["databases"],
                "bacteriome/taxdump.{assembler}.{dereper}"),
            db = os.path.join(config["output"]["databases"],
                "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index")
        conda:
            config["envs"]["krakenuniq"]
        threads:
            config["params"]["databases"]["threads"]
        shell:
            '''
            dbdir=$(realpath {params.db})
            rm -rf $dbdir
            mkdir -p $dbdir

            fnadir=$(realpath $(dirname {input.fnadone}))
            dmpdir=$(realpath {params.tax})

            cp {input.taxdump[3]} $fnadir/

            pushd $dbdir
            ln -s $fnadir library
            ln -s $dmpdir taxonomy
            popd

            jellyfish=$(command -v jellyfish)
            export JELLYFISH_BIN=$jellyfish

            krakenuniq-build \
            --threads {threads} \
            --db $dbdir \
            >{log} 2>&1
            '''


    rule databases_bacteriome_krakenuniq_bracken_build:
        input:
            expand(os.path.join(config["output"]["databases"],
                "bacteriome/databases.{{assembler}}.{{dereper}}/krakenuniq/index/{ku}"),
                ku=["database.kdb", "database.idx", "taxDB"])
        output:
            expand(os.path.join(config["output"]["databases"],
                "bacteriome/databases.{{assembler}}.{{dereper}}/krakenuniq/index/database{kmers}mers.{suffix}"),
                kmers=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["kmers"],
                suffix=["kmer_distrib", "kraken"])
        log:
            os.path.join(config["output"]["databases"],
                "logs/krakenuniq/bracken_build.{assembler}.{dereper}.log")
        benchmark:
            os.path.join(config["output"]["databases"],
                "benchmark/krakenuniq/bracken_build.{assembler}.{dereper}.benchmark.txt")
        params:
            db = os.path.join(config["output"]["databases"],
                "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index"),
            ksize=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["ksize"],
            kmers=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["kmers"]
        conda:
            config["envs"]["krakenuniq"]
        threads:
            config["params"]["databases"]["threads"]
        shell:
            '''
            dbdir=$(realpath {params.db})

            rm -rf {log}

            for kmer in {params.kmers};
            do
                bracken-build -d $dbdir -t {threads} -k {params.ksize} -l $kmer >>{log} 2>&1
            done
            '''


    rule databases_bacteriome_krakenuniq_all:
        input:
            expand([
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index/{ku}"),
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index/database{kmers}mers.{suffix}")],
                    ku=["database.kdb", "database.idex", "taxDB"],
                    kmers=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["kmers"],
                    suffix=["kmer_distrib", "kraken"],
                    assembler=ASSEMBLERS,
                    dereper=DEREPERS)

else:
    rule databases_bacteriome_krakenuniq_all:
        input:


rule databases_bacteriome_all:
    input:
        rules.databases_bacteriome_kmcp_all.input,
        rules.databases_bacteriome_kraken2_all.input,
        rules.databases_bacteriome_krakenuniq_all.input


rule databases_all:
    input:
        rules.databases_bacteriome_all.input


localrules:
    databases_bacteriome_kmcp_all,
    databases_bacteriome_kraken2_all,
    databases_bacteriome_krakenuniq_all