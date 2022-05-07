if config["params"]["identify"]["virsorter2"]["do"]:
    rule identify_virsorter2_setup_db:
        output:
            os.path.join(config["params"]["identify"]["virsorter2"]["db"], "Done_all_setup")
        conda:
            config["envs"]["virsorter2"]
        benchmark:
            os.path.join(config["output"]["identify"], "benchmark/virsorter2_setup_db.txt")
        log:
            os.path.join(config["output"]["identify"], "logs/virsorter2_setup_db.log")
        threads:
            config["params"]["identify"]["threads"]
        params:
            db_dir = config["params"]["identify"]["virsorter2"]["db"]
        shell:
            '''
            mkdir -p {params.db_dir}

            virsorter setup --db-dir {params.db_dir} --jobs {threads} >{log} 2>&1
            '''


    rule identify_virsorter2_config:
        input:
            os.path.join(config["params"]["identify"]["virsorter2"]["db"], "Done_all_setup")
        output:
            os.path.join(config["params"]["identify"]["virsorter2"]["db"], "Done_all_config")
        log:
            os.path.join(config["output"]["identify"], "logs/virsorter2_config.log")
        conda:
            config["envs"]["virsorter2"]
        threads:
            config["params"]["identify"]["threads"]
        params:
            db_dir = config["params"]["identify"]["virsorter2"]["db"]
        shell:
            '''
            virsorter config --init-source --db-dir={params.db_dir} >{log} 2>&1

            virsorter config --set GENERAL_THREADS={threads} >>{log} 2>&1
            virsorter config --set HMMSEARCH_THREADS={threads} >>{log} 2>&1
            virsorter config --set CLASSIFY_THREADS={threads} >>{log} 2>&1

            touch {output}
            '''


    rule identify_virsorter2_run:
        input:
            config_done = os.path.join(config["params"]["identify"]["virsorter2"]["db"], "Done_all_config"),
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(config["output"]["identify"],
                                "virsorter2/{{assembly_group}}.{{assembler}}.vs2.out/final-viral-{suffix}"),
                   suffix=["combined.fa", "score.tsv", "boundary.tsv"])
        benchmark:
            os.path.join(config["output"]["identify"], "benchmark/virsorter2/virsorter2.{assembly_group}.{assembler}.benchmark.txt")
        log:
            os.path.join(config["output"]["identify"], "logs/virsorter2/virsorter2.{assembly_group}.{assembler}.log")
        conda:
            config["envs"]["virsorter2"]
        threads:
            config["params"]["identify"]["threads"]
        params:
            working_dir = os.path.join(config["output"]["identify"], "virsorter2/{assembly_group}.{assembler}.vs2.out"),
            include_groups = ",".join(config["params"]["identify"]["virsorter2"]["include_groups"]),
            min_length = config["params"]["identify"]["virsorter2"]["min_length"],
            min_score = config["params"]["identify"]["virsorter2"]["min_score"],
            provirus_off = "--provirus-off" if config["params"]["identify"]["virsorter2"]["provirus_off"] else "",
            prep_for_dramv = "--prep-for-dramv" if config["params"]["identify"]["virsorter2"]["prep_for_dramv"] else "",
            rm_tmpdir = "--rm-tmpdir" if config["params"]["identify"]["virsorter2"]["rm_tmpdir"] else ""
        shell:
            '''
            virsorter run \
            {params.prep_for_dramv} \
            {params.provirus_off} \
            {params.rm_tmpdir} \
            --working-dir {params.working_dir} \
            --seqfile {input.scaftigs} \
            --include-groups {params.include_groups} \
            --min-length {params.min_length} \
            --min-score {params.min_score} \
            --jobs {threads} all \
            >{log} 2>&1
            '''


    rule identify_virsorter2_run_all:
        input:
            expand(os.path.join(config["output"]["identify"],
                                "virsorter2/{assembly_group}.{assembler}.vs2.out/final-viral-{suffix}"),
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST,
                   assembler=ASSEMBLERS,
                   suffix=["combined.fa", "score.tsv", "boundary.tsv"])
    
else:
    rule identify_virsorter2_run_all:
        input:


rule identify_all:
    input:
        rules.identify_virsorter2_run_all.input


localrules:
    identify_virsorter2_setup_db,
    identify_virsorter2_config,
    identify_virsorter2_run_all,
    identify_all