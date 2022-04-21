if config["params"]["identify"]["virsorter2"]["do"]:
    rule identify_virsorter2_setup_db:
        output:
            os.path.join(config["params"]["identify"]["virsorter2"]["db"], "Done_all_setup")
        params:
            db_dir = config["params"]["identify"]["virsorter2"]["db"]
        conda:
            config["envs"]["virsorter2"]
        benchmark:
            os.path.join(config["output"]["identify"], "benchmark/virsorter2_setup_db.txt")
        log:
            os.path.join(config["output"]["identify"], "logs/virsorter2_setup_db.log")
        threads:
            config["params"]["identify"]["threads"]
        shell:
            '''
            mkdir -p {params.db_dir}

            virsorter setup -d {params.db_dir} -j {threads} >{log} 2>&1
            '''


    rule identify_virsorter2_run:
        input:

        output:

        shell:
            '''
            
            '''