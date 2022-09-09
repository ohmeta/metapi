rule predict_mags_gene_prodigal:
    input:
        binning_done = lambda wildcards: get_binning_done(wildcards, [wildcards.binner_checkm])
    output:
        predict_done = os.path.join(
            config["output"]["predict"],
            "mags_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binner_checkm}/predict_done")
    log:
        os.path.join(config["output"]["predict"],
                     "logs/mags_gene/{binning_group}.{assembly_group}.{assembler}.{binner_checkm}.prodigal.log")
    conda:
        config["envs"]["checkm"]
    params:
        wrapper_dir = WRAPPER_DIR
    threads:
        config["params"]["predict"]["threads"]
    shell:
        '''
        python {params.wrapper_dir}/prodigal_wrapper.py \
        {threads} \
        {input.binning_done} \
        {output.predict_done}
        '''


rule predict_mags_gene_prodigal_report:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "mags_gene/{binning_group}.{assembly_group}.{{assembler}}.prodigal/{{binner_checkm}}/predict_done"),
                zip,
                binning_group=ASSEMBLY_GROUP["binning_group"],
                assembly_group=ASSEMBLY_GROUP["assembly_group"])
        output:
            os.path.join(config["output"]["predict"],
                         "report/mags_gene_stats_{assembler}_{binner_checkm}.tsv")
        run:
            import pandas as pd

            shell(f'''mkdir -p {os.path.dirname(output[0])}''')

            df_list = [pd.read_csv(i, sep="\t") for i in input]
            pd.concat(df_list, axis=0).to_csv(output[0], sep="\t", index=False)

 
rule predict_mags_gene_prodigal_all:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "report/mags_gene_stats_{assembler}_{binner_checkm}.tsv"),
            assembler=ASSEMBLERS,
            binner_checkm=BINNERS_CHECKM)

        #rules.binning_all.input


if config["params"]["predict"]["mags_to_gene"]["prokka"]["do"]:
    rule predict_mags_gene_prokka:
        input:
            mags_dir = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/{binner_checkm}")
        output:
            done = os.path.join(
                config["output"]["predict"],
                "mags_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{binner_checkm}/predict_done")
        conda:
            config["envs"]["predict"]
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "mags_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{binner_checkm}"),
            kingdom = config["params"]["predict"]["mags_to_gene"]["prokka"]["kingdom"]
        log:
            os.path.join(config["output"]["predict"],
                         "logs/mags_gene/{binning_group}.{assembly_group}.{assembler}.{binner_checkm}.prokka.log")
        threads:
            config["params"]["predict"]["threads"]
        script:
            "../wrappers/prokka_wrapper.py"


    rule predict_mags_gene_prokka_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["predict"],
                    "mags_gene/{binning_group}.{assembly_group}.{{assembler}}.prokka/{{binner_checkm}}/predict_done"),
                    zip,
                    binning_group=ASSEMBLY_GROUP["binning_group"],
                    assembly_group=ASSEMBLY_GROUP["assembly_group"])
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/mags_gene_{assembler}.{binner_checkm}.multiqc/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/mags_gene_{assembler}.{binner_checkm}.multiqc/prokka_multiqc_report_data"))
        conda:
            config["envs"]["multiqc"]
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/report/mags_gene_{assembler}.{binner_checkm}.multiqc.prokka.log")
        params:
            input_dir = lambda wildcards: expand(expand(os.path.join(
                config["output"]["predict"],
                "mags_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{{binner_checkm}}"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
                binner_checkm=wildcards.binner_checkm),
            output_dir = os.path.join(
                config["output"]["predict"],
                "report/mags_gene_{assembler}.{binner_checkm}.multiqc")
        shell:
            '''
            multiqc \
            --cl_config "prokka_fn_snames: True" \
            --outdir {params.output_dir} \
            --title prokka \
            --module prokka \
            {params.input_dir} \
            2> {log}
            '''


    rule predict_mags_gene_prokka_all:
        input:
            expand(expand(
                os.path.join(
                    config["output"]["predict"],
                    "mags_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{binner_checkm}/predict_done"),
                    zip,
                    binning_group=ASSEMBLY_GROUPS["binning_group"],
                    assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                    assembler=ASSEMBLY_GROUPS["assembler"]),
                    binner_checkm=BINNERS_CHECKM),
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "report/mags_gene_{assembler}.{binner_checkm}.multiqc/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/mags_gene_{assembler}.{binner_checkm}.multiqc/prokka_multiqc_report_data")],
                   assembler=ASSEMBLERS,
                   binner_checkm=BINNERS_CHECKM)

            #rules.binning_all.input

else:
    rule predict_mags_gene_prokka_all:
        input:


rule predict_mags_gene_all:
    input:
        rules.predict_mags_gene_prodigal_all.input,
        rules.predict_mags_gene_prokka_all.input,


rule predict_all:
    input:
        rules.predict_scaftigs_gene_all.input,
        rules.predict_mags_gene_all.input


localrules:
    predict_mags_gene_prodigal_report,
    predict_mags_gene_prokka_all,
    predict_mags_gene_all,
    predict_all