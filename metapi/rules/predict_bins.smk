rule predict_bins_gene_prodigal:
    input:
        binning_done = lambda wildcards: get_binning_done(wildcards, [wildcards.binner_checkm])
    output:
        predict_done = os.path.join(
            config["output"]["predict"],
            "bins_gene/{assembly_group}.{assembler}.prodigal.out/{binner_checkm}/predict_done")
    log:
        os.path.join(config["output"]["predict"],
                     "logs/bins_gene/{assembly_group}.{assembler}.{binner_checkm}.prodigal.log")
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


rule predict_bins_gene_prodigal_report:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembly_group}.{{assembler}}.prodigal.out/{{binner_checkm}}/predict_done"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
        output:
            os.path.join(config["output"]["predict"],
                         "report/bins_gene_stats_{assembler}_{binner_checkm}.tsv")
        run:
            import pandas as pd

            shell(f'''mkdir -p {os.path.dirname(output[0])}''')

            df_list = [pd.read_csv(i, sep="\t") for i in input]
            pd.concat(df_list, axis=0).to_csv(output[0], sep="\t", index=False)

 
rule predict_bins_gene_prodigal_all:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "report/bins_gene_stats_{assembler}_{binner_checkm}.tsv"),
            assembler=ASSEMBLERS,
            binner_checkm=BINNERS_CHECKM)

        #rules.binning_all.input


if config["params"]["predict"]["bins_to_gene"]["prokka"]["do"]:
    rule predict_bins_gene_prokka:
        input:
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/{binner_checkm}")
        output:
            done = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembly_group}.{assembler}.prokka.out/{binner_checkm}/predict_done")
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembly_group}.{assembler}.prokka.out/{binner_checkm}"),
            kingdom = config["params"]["predict"]["bins_to_gene"]["prokka"]["kingdom"]
        log:
            os.path.join(config["output"]["predict"],
                         "logs/bins_gene/{assembly_group}.{assembler}.{binner_checkm}.prokka.log")
        threads:
            config["params"]["predict"]["threads"]
        run:
            import glob
            import os
            import time
            import subprocess

            bin_list = glob.glob(input.bins_dir + "/*bin*fa")
            gff_count = 0

            for bin_fa in bin_list:
                bin_id = os.path.basename(os.path.splitext(bin_fa)[0])
                output_dir = os.path.join(params.output_dir, bin_id)
                gff_file = os.path.join(output_dir, bin_id + ".gff")

                shell(
                    f'''
                    echo "\nProcessing {bin_fa}\n" >> {log}

                    prokka {bin_fa} \
                    --force \
                    --centre X \
                    --compliant \
                    --cpus {threads} \
                    --outdir {output_dir} \
                    --locustag {bin_id} \
                    --prefix {bin_id} \
                    --kingdom {params.kingdom} \
                    2>> {log} 
                    ''')

                if os.path.exists(gff_file):
                    gff_count += 1

            if gff_count == len(bin_list):
                shell('''touch {output.done}''')


    rule predict_bins_gene_prokka_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene/{assembly_group}.{{assembler}}.prokka.out/{{binner_checkm}}/predict_done"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/bins_gene_{assembler}.{binner_checkm}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/bins_gene_{assembler}.{binner_checkm}.multiqc.out/prokka_multiqc_report_data"))
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/report/bins_gene_{assembler}.{binner_checkm}.multiqc.prokka.log")
        params:
            input_dir = lambda wildcards: expand(os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembly_group}.{assembler}.prokka.out/{binner_checkm}"),
                assembler=wildcards.assembler, 
                binner_checkm=wildcards.binner_checkm,
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST),
            output_dir = os.path.join(
                config["output"]["predict"],
                "report/bins_gene_{assembler}.{binner_checkm}.multiqc.out")
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


    rule predict_bins_gene_prokka_all:
        input:
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene/{assembly_group}.{assembler}.prokka.out/{binner_checkm}/predict_done"),
                os.path.join(
                    config["output"]["predict"],
                    "report/bins_gene_{assembler}.{binner_checkm}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/bins_gene_{assembler}.{binner_checkm}.multiqc.out/prokka_multiqc_report_data")],
                   assembler=ASSEMBLERS,
                   binner_checkm=BINNERS_CHECKM,
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)#,

            #rules.binning_all.input

else:
    rule predict_bins_gene_prokka_all:
        input:


rule predict_bins_gene_all:
    input:
        rules.predict_bins_gene_prodigal_all.input,
        rules.predict_bins_gene_prokka_all.input,


rule predict_all:
    input:
        rules.predict_scaftigs_gene_all.input,
        rules.predict_bins_gene_all.input


localrules:
    predict_bins_gene_prodigal_report,
    predict_bins_gene_prokka_all,
    predict_bins_gene_all,
    predict_all