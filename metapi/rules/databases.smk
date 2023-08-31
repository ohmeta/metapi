rule databases_bacteriome_refine_taxonomy:
    input:
        rep_genomes_info = os.path.join(
            config["output"]["dereplicate"],
            "report/bacteriome/checkm_table_genomes_info.{assembler}.{dereper}.tsv.gz"),
        table_gtdb = os.path.join(
            config["output"]["taxonomic"],
            "report/gtdbtk/MAGs_hmq_{assembler}_all_gtdbtk_gtdb.tsv")
    output:
        taxonomy = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxonomy.{assembler}.{dereper}/taxonomy.tsv"),
        scaftigs = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxonomy.{assembler}.{dereper}/scaftigs.tsv"),
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
            output.taxonomy,
            output.scaftigs)

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
            taxdump=["names.dmp", "nodes.dmp", "taxid.map"]),
        taxmap = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/data/{taxmap}.map"),
            taxmap=["taxid", "name"])
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_generate_taxdump/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_generate_taxdump/{assembler}.{dereper}.txt")
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

        mkdir -p {params.out_dir}/data

        cp {params.out_dir}/taxid.map {params.out_dir}/data/taxid.map

        cat {params.out_dir}/taxid.map | \
        taxonkit lineage \
        --data-dir {params.out_dir} \
        -i 2 -n -L | \
        cut -f 1,3 > {params.out_dir}/data/name.map 2>> {log}
        '''


rule databases_bacteriome_generate_prelim_map:
    input:
        fnadone = os.path.join(config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/done"),
        scaftigs = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxonomy.{assembler}.{dereper}/scaftigs.tsv"),
        taxid = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/taxid.map")
    output:
        prelim_map = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/prelim_map.txt"),
        prelim_map_fna = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/prelim_map.txt")
    run:
        import os
        import pandas as pd

        scaftigs_info = pd.read_csv(input.scaftigs, sep="\t", names=["genome_id", "scaftigs_id"])
        taxid_info = pd.read_csv(input.taxid, sep="\t", names=["genome_id", "tax_id"])
        prelim_map = pd.merge(scaftigs_info, taxid_info, how="inner", on=["genome_id"])
        prelim_map["first"] = "TAXID"

        prelim_map.loc[:, ["first", "scaftigs_id", "tax_id"]]\
        .to_csv(output.prelim_map, sep="\t", index=False, header=None)

        fnadir = os.path.dirname(input.fnadone)
        shell("cp {output.prelim_map} {output.prelim_map_fna}")


rule databases_bacteriome_extract_taxonomy:
    input:
        taxdump = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
            taxdump=["names.dmp", "nodes.dmp"])
    output:
        taxtab = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab")
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_extract_taxonomy/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_extract_taxonomy/{assembler}.{dereper}.txt")
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
        fnadone = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/done")
    output:
        directory(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kmcp/compute"))
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_kmcp_compute/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_kmcp_compute/{assembler}.{dereper}.txt")
    params:
        kmer = config["params"]["databases"]["bacteriome"]["kmcp"]["compute"]["kmer"],
        split_number = config["params"]["databases"]["bacteriome"]["kmcp"]["compute"]["split_number"]
    threads:
        config["params"]["databases"]["threads"]
    conda:
        config["envs"]["kmcp"]
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
        taxmap = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/data/{taxmap}.map"),
            taxmap=["taxid", "name"]),
        compute = os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kmcp/compute")
    output:
        directory(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kmcp/index"))
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_kmcp_index/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_kmcp_index/{assembler}.{dereper}.txt")
    params:
        false_positive_rate = config["params"]["databases"]["bacteriome"]["kmcp"]["index"]["false_positive_rate"]
    threads:
        config["params"]["databases"]["threads"]
    conda:
        config["envs"]["kmcp"]
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


if config["params"]["databases"]["bacteriome"]["kmcp"]["do"]:
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


rule databases_bacteriome_kraken2_build:
    input:
        fnadone = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/done"),
        taxtab = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab"),
        taxdump = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
            taxdump=["names.dmp", "nodes.dmp"]),
        prelim_map = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/prelim_map.txt"),
        prelim_map_fna = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/prelim_map.txt")
    output:
        expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{{assembler}}.{{dereper}}/kraken2/index/{k2}.k2d"),
            k2=["hash", "opts", "taxo"])
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_kraken2_build/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_kraken2_build/{assembler}.{dereper}.txt")
    params:
        tax = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}"),
        db = os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kraken2/index")
    threads:
        config["params"]["databases"]["threads"]
    conda:
        config["envs"]["kraken2"]
    shell:
        '''
        dbdir=$(realpath {params.db})
        rm -rf $dbdir
        mkdir -p $dbdir

        fnadir=$(realpath $(dirname {input.fnadone}))
        dmpdir=$(realpath {params.tax})

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
        expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{{assembler}}.{{dereper}}/kraken2/index/{k2}.k2d"),
            k2=["hash", "opts", "taxo"])
    output:
        expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{{assembler}}.{{dereper}}/kraken2/index/database{kmers}mers.{suffix}"),
            kmers=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["kmers"],
            suffix=["kmer_distrib", "kraken"])
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_kraken2_bracken_build/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_kraken2_bracken_build/{assembler}.{dereper}.txt")
    params:
        db = os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/kraken2/index"),
        ksize=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["ksize"],
        kmers=config["params"]["databases"]["bacteriome"]["kraken2"]["bracken"]["kmers"]
    threads:
        config["params"]["databases"]["threads"]
    conda:
        config["envs"]["kraken2"]
    shell:
        '''
        dbdir=$(realpath {params.db})

        rm -rf {log}

        for kmer in {params.kmers};
        do
            bracken-build -d $dbdir -t {threads} -k {params.ksize} -l $kmer >>{log} 2>&1
        done
        '''


if config["params"]["databases"]["bacteriome"]["kraken2"]["do"]:
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


rule databases_bacteriome_krakenuniq_build:
    input:
        fnadone = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/done"),
        taxtab = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/taxonomy.tab"),
        taxdump = expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{{assembler}}.{{dereper}}/{taxdump}"),
            taxdump=["names.dmp", "nodes.dmp"]),
        prelim_map = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}/prelim_map.txt"),
        prelim_map_fna = os.path.join(
            config["output"]["databases"],
            "bacteriome/fna.{assembler}.{dereper}/prelim_map.txt")
    output:
        expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{{assembler}}.{{dereper}}/krakenuniq/index/{ku}"),
            ku=["database.kdb", "database.idx", "taxDB"])
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_krakenuniq_build/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_krakenuniq_build/{assembler}.{dereper}.txt")
    params:
        tax = os.path.join(
            config["output"]["databases"],
            "bacteriome/taxdump.{assembler}.{dereper}"),
        db = os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index")
    threads:
        config["params"]["databases"]["threads"]
    conda:
        config["envs"]["krakenuniq"]
    shell:
        '''
        dbdir=$(realpath {params.db})
        rm -rf $dbdir
        mkdir -p $dbdir

        fnadir=$(realpath $(dirname {input.fnadone}))
        dmpdir=$(realpath {params.tax})

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
        expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{{assembler}}.{{dereper}}/krakenuniq/index/{ku}"),
            ku=["database.kdb", "database.idx", "taxDB"])
    output:
        expand(os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{{assembler}}.{{dereper}}/krakenuniq/index/database{kmers}mers.{suffix}"),
            kmers=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["kmers"],
            suffix=["kmer_distrib", "kraken"])
    log:
        os.path.join(
            config["output"]["databases"],
            "logs/databases_bacteriome_krakenuniq_bracken_build/{assembler}.{dereper}.log")
    benchmark:
        os.path.join(
            config["output"]["databases"],
            "benchmark/databases_bacteriome_krakenuniq_bracken_build/{assembler}.{dereper}.txt")
    params:
        db = os.path.join(
            config["output"]["databases"],
            "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index"),
        ksize=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["ksize"],
        kmers=config["params"]["databases"]["bacteriome"]["krakenuniq"]["bracken"]["kmers"]
    threads:
        config["params"]["databases"]["threads"]
    conda:
        config["envs"]["krakenuniq"]
    shell:
        '''
        dbdir=$(realpath {params.db})

        rm -rf {log}

        for kmer in {params.kmers};
        do
            bracken-build -d $dbdir -t {threads} -k {params.ksize} -l $kmer >>{log} 2>&1
        done
        '''


if config["params"]["databases"]["bacteriome"]["krakenuniq"]["do"]:
    rule databases_bacteriome_krakenuniq_all:
        input:
            expand([
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index/{ku}"),
                os.path.join(
                    config["output"]["databases"],
                    "bacteriome/databases.{assembler}.{dereper}/krakenuniq/index/database{kmers}mers.{suffix}")],
                    ku=["database.kdb", "database.idx", "taxDB"],
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