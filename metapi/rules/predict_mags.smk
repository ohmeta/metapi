rule predict_mags_gene_prodigal:
    input:
        binning_done = lambda wildcards: get_binning_done(wildcards, [wildcards.binner_checkm])
    output:
        predict_done = os.path.join(
            config["output"]["predict"],
            "mags_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binner_checkm}/predict_done")
    log:
        os.path.join(
            config["output"]["predict"],
            "logs/predict_mags_gene_prodigal/{binning_group}.{assembly_group}.{assembler}.{binner_checkm}.log")
    benchmark:
        os.path.join(
            config["output"]["predict"],
            "benchmark/predict_mags_gene_prodigal/{binning_group}.{assembly_group}.{assembler}.{binner_checkm}.txt")
    params:
        wrapper_dir = WRAPPER_DIR
    threads:
        config["params"]["predict"]["threads"]
    conda:
        config["envs"]["checkm"]
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
            os.path.join(config["output"]["predict"], "report/mags_gene_stats_{assembler}_{binner_checkm}.tsv.gz")
        run:
            import pandas as pd

            shell(f'''mkdir -p {os.path.dirname(output[0])}''')

            df_list = [pd.read_csv(i, sep="\t") for i in input]
            pd.concat(df_list, axis=0).to_csv(output[0], sep="\t", index=False)


checkm_df_list = []
for binner_checkm in BINNERS_CHECKM:
    checkm_df = ASSEMBLY_GROUPS.copy()
    checkm_df["binner_checkm"] = binner_checkm 
    checkm_df_list.append(checkm_df)
CHECKM_GROUPS = pd.concat(checkm_df_list, axis=0)


rule predict_mags_gene_prodigal_all:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "mags_gene/{binning_group}.{assembly_group}.{assembler}.prodigal/{binner_checkm}/predict_done"),
            zip,
            binning_group=CHECKM_GROUPS["binning_group"],
            assembly_group=CHECKM_GROUPS["assembly_group"],
            assembler=CHECKM_GROUPS["assembler"],
            binner_checkm=CHECKM_GROUPS["binner_checkm"]),
        expand(os.path.join(
            config["output"]["predict"],
            "report/mags_gene_stats_{assembler}_{binner_checkm}.tsv.gz"),
            assembler=ASSEMBLERS,
            binner_checkm=BINNERS_CHECKM)


rule predict_mags_gene_prokka:
    input:
        mags_dir = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/{binner_checkm}")
    output:
        done = os.path.join(
            config["output"]["predict"],
            "mags_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{binner_checkm}/predict_done")
    log:
        os.path.join(
            config["output"]["predict"],
            "logs/predict_mags_gene_prokka/{binning_group}.{assembly_group}.{assembler}.{binner_checkm}.log")
    benchmark:
        os.path.join(
            config["output"]["predict"],
            "benchmark/predict_mags_gene_prokka/{binning_group}.{assembly_group}.{assembler}.{binner_checkm}.txt")
    params:
        output_dir = os.path.join(
            config["output"]["predict"],
            "mags_gene/{binning_group}.{assembly_group}.{assembler}.prokka/{binner_checkm}"),
        kingdom = config["params"]["predict"]["mags_to_gene"]["prokka"]["kingdom"]
    threads:
        config["params"]["predict"]["threads"]
    conda:
        config["envs"]["predict"]
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
    log:
        os.path.join(
            config["output"]["predict"],
            "logs/predict_mags_gene_prokka_multiqc/{assembler}.{binner_checkm}.log")
    benchmark:
        os.path.join(
            config["output"]["predict"],
            "benchmark/predict_mags_gene_prokka_multiqc/{assembler}.{binner_checkm}.txt")
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
    conda:
        config["envs"]["multiqc"]
    shell:
        '''
        multiqc \
        --cl_config "prokka_fn_snames: True" \
        --outdir {params.output_dir} \
        --title prokka \
        --module prokka \
        {params.input_dir} \
        >{log} 2>&1
        '''


if config["params"]["predict"]["mags_to_gene"]["prokka"]["do"]:
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

else:
    rule predict_mags_gene_prokka_all:
        input:


rule predict_mags_gene_all:
    input:
        rules.predict_mags_gene_prodigal_all.input,
        rules.predict_mags_gene_prokka_all.input


rule predict_all:
    input:
        rules.predict_scaftigs_gene_all.input,
        rules.predict_mags_gene_all.input


localrules:
    predict_mags_gene_prodigal_report,
    predict_mags_gene_prokka_all,
    predict_mags_gene_all,
    predict_all