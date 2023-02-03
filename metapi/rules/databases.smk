rule databases_bacteriome_refine_taxonomy:
    input:
        genomes_info = os.path.join(
            config["output"]["check"],
            "report/checkm/checkm_table_genomes_info.{assembler}.all.tsv"),
        table_gtdb = os.path.join(
            config["output"]["taxonomic"],
            "report/gtdbtk/MAGs_hmq.rep.{assembler}.{dereper}.gtdbtk.gtdb.tsv")
    output:
        taxonomy = os.path.join(config["output"]["databases"], "bacteriome/taxonomy.{assembler}.{dereper}/taxonomy.tsv")
    params:
        db_name = config["params"]["databases"]["bacteriome"]["name"],
        rep_level = config["params"]["databases"]["bacteriome"]["rep_level"],
        out_dir = os.path.join(config["output"]["databases"], "bacteriome"),
        base_dir = os.path.join(config["output"]["databases"], "bacteriome/fna.{assembler}.{dereper}")
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


rule databases_all:
    input:
        expand(
            os.path.join(config["output"]["databases"], "bacteriome/taxonomy.{assembler}.{dereper}/taxonomy.tsv"),
            assembler=ASSEMBLERS,
            dereper=DEREPERS)