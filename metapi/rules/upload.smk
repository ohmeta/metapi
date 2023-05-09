if config["upload"]["do"]:
    rule upload_generate_samples_info:
        input:
            config["params"]["samples"]
        output:
            os.path.join(config["output"]["upload"], "table/MIxS_Samples.xlsx")
        run:
            metapi.gen_samples_info(SAMPLES, output[0], config)


    rule upload_md5_short_reads:
        input:
            alignment_input_with_short_reads
        output:
            os.path.join(config["output"]["upload"], "md5/short_reads/{sample}.md5")
        shell:
            '''
            md5sum {input} > {output}
            '''


    rule upload_generate_run_info:
        input:
            expand(os.path.join(
                config["output"]["upload"], "md5/short_reads/{sample}.md5"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["upload"], "table/Experiment_Run.xlsx")
        threads:
            config["upload"]["threads"]
        run:
            metapi.gen_info(input, output[0], config, threads, "sequencing_run")


    rule upload_sequencing_all:
        input:
            os.path.join(config["output"]["upload"], "table/Experiment_Run.xlsx"),
            os.path.join(config["output"]["upload"], "table/MIxS_Samples.xlsx")


    localrules:
        upload_generate_samples_info,
        upload_generate_run_info,


    if len(ASSEMBLERS) != 0:
        rule upload_md5_scaftigs:
            input:
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
            output:
                os.path.join(
                    config["output"]["upload"],
                    "md5/scaftigs/{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.md5")
            shell:
                '''
                md5sum {input} > {output}
                '''


        rule upload_generate_assembly_info:
            input:
                expand(os.path.join(
                    config["output"]["upload"],
                    "md5/scaftigs/{{assembler}}/{binning_group}.{assembly_group}.{{assembler}}.scaftigs.md5"),
                    zip,
                    binning_group=ASSEMBLY_GROUP["binning_group"],
                    assembly_group=ASSEMBLY_GROUP["assembly_group"])
            output:
                os.path.join(
                    config["output"]["upload"],
                    "table/Genome_Assembly_{assembler}.xlsx")
            threads:
                config["upload"]["threads"]
            run:
                metapi.gen_info(input, output[0], config, threads, "assembly")


        rule upload_assembly_all:
            input:
                expand(os.path.join(
                    config["output"]["upload"],
                    "table/Genome_Assembly_{assembler}.xlsx"),
                    assembler=ASSEMBLERS)


        localrules:
            upload_generate_assembly_info


    else:
        rule upload_assembly_all:
            input:

else:
    rule upload_sequencing_all:
        input:


    rule upload_assembly_all:
        input:


rule upload_all:
    input:
        rules.upload_sequencing_all.input,
        rules.upload_assembly_all.input#,


localrules:
    upload_sequencing_all,
    upload_assembly_all,
    upload_all