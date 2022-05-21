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

# https://github.com/EddyRivasLab/hmmer/issues/161
# hmmsearch threads: 2 (recommand)
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


    rule identify_virsorter2:
        input:
            config_done = os.path.join(config["params"]["identify"]["virsorter2"]["db"], "Done_all_config"),
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(config["output"]["identify"],
                                "virsorter2/{{assembly_group}}.{{assembler}}.vs2.out/{{assembly_group}}.{{assembler}}-final-viral-{suffix}"),
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
            label = "{assembly_group}.{assembler}",
            working_dir = os.path.join(config["output"]["identify"], "virsorter2/{assembly_group}.{assembler}.vs2.out"),
            include_groups = ",".join(config["params"]["identify"]["virsorter2"]["include_groups"]),
            min_length = config["params"]["identify"]["virsorter2"]["min_length"],
            min_score = config["params"]["identify"]["virsorter2"]["min_score"],
            provirus_off = "--provirus-off" if config["params"]["identify"]["virsorter2"]["provirus_off"] else "",
            prep_for_dramv = "--prep-for-dramv" if config["params"]["identify"]["virsorter2"]["prep_for_dramv"] else "",
            rm_tmpdir = "--rm-tmpdir" if config["params"]["identify"]["virsorter2"]["rm_tmpdir"] else "",
            keep_original_seq = "--keep-original-seq" if config["params"]["identify"]["virsorter2"]["keep_original_seq"] else ""
        shell:
            '''
            set +e

            virsorter run \
            {params.prep_for_dramv} \
            {params.provirus_off} \
            {params.rm_tmpdir} \
            {params.keep_original_seq} \
            --working-dir {params.working_dir} \
            --seqfile {input.scaftigs} \
            --label {params.label} \
            --include-groups {params.include_groups} \
            --min-length {params.min_length} \
            --min-score {params.min_score} \
            --jobs {threads} all \
            >{log} 2>&1

            exitcode=$?
            echo "Exit code is: $exitcode" >> {log}

            if [ $exitcode -eq 1 ];
            then
                grep -oEi "No genes from the contigs are left in iter-0/all.pdg.faa after preprocess" {log} 
                grepcode=$?
                if [ $grepcode -eq 0 ];
                then
                    touch {output[0]} >> {log} 2>&1
                    touch {output[1]} >> {log} 2>&1
                    touch {output[2]} >> {log} 2>&1
                    echo "Touch empty file: {output[0]}" >> {log} 2>&1
                    echo "Touch empty file: {output[1]}" >> {log} 2>&1
                    echo "Touch empty file: {output[2]}" >> {log} 2>&1
                    exit 0
                else
                    grep -oEi "Error in rule circular_linear_split" {log}
                    grepcode=$?
                    if [ $grepcode -eq 0 ];
                    then
                        touch {output[0]} >> {log} 2>&1
                        touch {output[1]} >> {log} 2>&1
                        touch {output[2]} >> {log} 2>&1
                        echo "Touch empty file: {output[0]}" >> {log} 2>&1
                        echo "Touch empty file: {output[1]}" >> {log} 2>&1
                        echo "Touch empty file: {output[2]}" >> {log} 2>&1
                        exit 0
                    else
                        echo "Runing failed, check Virsorter log please." >> {log} 2>&1
                        exit $exitcode
                    fi
                fi
            fi
            '''


    rule identify_virsorter2_all:
        input:
            expand(os.path.join(config["output"]["identify"],
                                "virsorter2/{assembly_group}.{assembler}.vs2.out/{assembly_group}.{assembler}-final-viral-{suffix}"),
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST,
                   assembler=ASSEMBLERS,
                   suffix=["combined.fa", "score.tsv", "boundary.tsv"])
 
else:
    rule identify_virsorter2_all:
        input:


if config["params"]["identify"]["deepvirfinder"]["do"]:
    rule identify_deepvirfinder:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            expand(
                os.path.join(
                    config["output"]["identify"],
                    "deepvirfinder/{{assembly_group}}.{{assembler}}.dvf.out/{{assembly_group}}.{{assembler}}.scaftigs.fa.gz_gt{min_length}bp_dvfpred.txt"),
                min_length=config["params"]["identify"]["deepvirfinder"]["min_length"])
        benchmark:
            os.path.join(config["output"]["identify"], "benchmark/deepvirfinder/deepvirfinder.{assembly_group}.{assembler}.benchmark.txt")
        log:
            os.path.join(config["output"]["identify"], "logs/deepvirfinder/deepvirfinder.{assembly_group}.{assembler}.log")
        conda:
            config["envs"]["deepvirfinder"]
        params:
            deepvirfinder = config["params"]["identify"]["deepvirfinder"]["script"],
            out_dir = os.path.join(config["output"]["identify"], "deepvirfinder/{assembly_group}.{assembler}.dvf.out"),
            min_length=config["params"]["identify"]["deepvirfinder"]["min_length"]
        threads:
            config["params"]["identify"]["threads"]
        shell:
            '''
            set +e

            python {params.deepvirfinder} \
            --in {input} \
            --out {params.out_dir} \
            --len {params.min_length} \
            --core {threads} \
            >{log} 2>&1

            exitcode=$?
            echo "Exit code is: $exitcode" >> {log} 2>&1

            if [ $exitcode -eq 1 ];
            then
                grep -oEi "ValueError: not enough values to unpack" {log} 
                grepcode=$?
                if [ $grepcode -eq 0 ];
                then
                    touch {output} >> {log} 2>&1
                    echo "Touch empty file: {output}" >> {log} 2>&1
                    exit 0
                else
                    echo "Runing failed, check DeepVirFinder log please." >> {log} 2>&1
                    exit $exitcode
                fi
            fi
            '''
 
    
    rule identify_deepvirfinder_all:
        input:
            expand(
                os.path.join(
                    config["output"]["identify"],
                    "deepvirfinder/{assembly_group}.{assembler}.dvf.out/{assembly_group}.{assembler}.scaftigs.fa.gz_gt{min_length}bp_dvfpred.txt"),
                min_length=config["params"]["identify"]["deepvirfinder"]["min_length"],
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST,
                assembler=ASSEMBLERS)

else:
    rule identify_deepvirfinder_all:
        input:


rule identify_single_all:
    input:
        rules.identify_virsorter2_all.input,
        rules.identify_deepvirfinder_all.input


localrules:
    identify_virsorter2_setup_db,
    identify_virsorter2_config,
    identify_virsorter2_all,
    identify_deepvirfinder_all,
    identify_single_all