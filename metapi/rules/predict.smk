PROKKA_SUFFIX = ["err", "log", "faa", "ffn", "fna", "fsa",
                 "gbk", "gff", "sqn", "tbl", "tsv", "txt"]


if config["params"]["predict"]["scaftigs_to_gene"]["prodigal"]["do"]:
    rule predict_scaftigs_gene_prodigal:
        input:
            os.path.join(config["output"]["assembly"],
                         "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            pep = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.faa"),
            cds = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.ffn"),
            gff = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.gff")
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{sample}.{assembler}.prodigal.log")
        params:
            format = config["params"]["predict"]["format"]
        shell:
            '''
            zcat {input} | \
            prodigal \
            -m \
            -a {output.pep} \
            -d {output.cds} \
            -o {output.gff} \
            -f {params.format} \
            -p meta -q \
            2> {log}
            '''


    rule predict_scaftigs_gene_prodigal_all:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.{ext}"),
                   ext=["faa", "ffn", "gff"],
                   assembler=ASSEMBLERS,
                   sample=SAMPLES.index.unique())

else:
    rule predict_scaftigs_gene_prodigal_all:
        input:


if config["params"]["predict"]["scaftigs_to_gene"]["prokka"]["do"]:
    rule predict_scaftigs_gene_prokka:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            expand(os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{{sample}}.{{assembler}}.prokka.out/{{sample}}.{{assembler}}.{ext}"),
                   ext=PROKKA_SUFFIX)
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{sample}.{assembler}.prokka.log")
        params:
            prefix = "{sample}.{assembler}",
            output_dir = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prokka.out")
        threads:
            config["params"]["predict"]["threads"]
        shell:
            '''
            prokka {input} \
            --force \
            --centre X \
            --compliant \
            --cpus {threads} \
            --outdir {params.output_dir} \
            --locustag {params.prefix} \
            --prefix {params.prefix} \
            --metagenome \
            2> {log}
            '''


    rule predict_scaftigs_gene_prokka_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{{assembler}}.prokka.out/{sample}.{{assembler}}.{ext}"),
                ext=PROKKA_SUFFIX,
                sample=SAMPLES.index.unique())
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data"))
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/report/scaftigs_gene_{assembler}.multiqc.prokka.log")
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "report/scaftigs_gene_{assembler}.multiqc.out")
        shell:
            '''
            multiqc \
            --cl_config "prokka_fn_snames: True" \
            --outdir {params.output_dir} \
            --title prokka \
            --module prokka \
            {input} \
            2> {log}
            '''


    rule predict_scaftigs_gene_prokka_all:
        input:
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{assembler}.prokka.out/{sample}.{assembler}.{ext}"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/scaftigs_gene_{assembler}.multiqc.out/prokka_multiqc_report_data")],
                   ext=PROKKA_SUFFIX,
                   assembler=ASSEMBLERS,
                   sample=SAMPLES.index.unique())

else:
    rule predict_scaftigs_gene_prokka_all:
        input:


if config["params"]["predict"]["bins_to_gene"]["prodigal"]["do"]:
    rule predict_bins_gene_prodigal:
        input:
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{sample}.{assembler}.out/{binner}")
        output:
            done = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prodigal.out/{sample}/done")
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prodigal.out/{sample}"),
            logs_dir = os.path.join(config["output"]["predict"],
                                    "logs/bins_gene/{sample}.prodigal"),
            format = config["params"]["predict"]["format"]
        run:
            import glob
            import os
            import time
            import subprocess

            bin_list = glob.glob(input.bins_dir + "/*bin*fa")
            gff_count = 0

            os.makedirs(params.output_dir, exist_ok=True)
            os.makedirs(params.logs_dir, exist_ok=True)

            for bin_fa in bin_list:
                bin_id = os.path.basename(os.path.splitext(bin_fa)[0])
                pep_file = os.path.join(params.output_dir, bin_id + ".faa")
                cds_file = os.path.join(params.output_dir, bin_id + ".ffn")
                gff_file = os.path.join(params.output_dir, bin_id + ".gff")
                log_file= os.path.join(params.logs_dir, bin_id + ".prodigal.log")

                shell(
                    '''
                    prodigal \
                    -i %s \
                    -m \
                    -a %s \
                    -d %s \
                    -o %s \
                    -f {params.format} \
                    -p single \
                    2> %s
                    ''' % (bin_fa, pep_file, cds_file, gff_file, log_file))

                if os.path.exists(gff_file):
                    gff_count += 1

            if gff_count == len(bin_list):
                shell('''touch {output.done}''')

       
    rule predict_bins_gene_prodigal_all:
        input:
            expand(os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prodigal.out/{sample}/done"),
                   assembler=ASSEMBLERS,
                   binner=BINNERS,
                   sample=SAMPLES.index.unique())

else:
    rule predict_bins_gene_prodigal_all:
        input:


if config["params"]["predict"]["bins_to_gene"]["prokka"]["do"]:
    rule predict_bins_gene_prokka:
        input:
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{sample}.{assembler}.out/{binner}")
        output:
            done = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prokka.out/{sample}/done")
        params:
            output_dir = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prokka.out/{sample}"),
            logs_dir = os.path.join(config["output"]["predict"],
                                    "logs/bins_gene/{sample}.prodigal"),
            kingdom = config["params"]["predict"]["bins_to_gene"]["prokka"]["kingdom"]
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
                log_file= os.path.join(params.logs_dir, bin_id + ".prokka.log")
                os.makedirs(params.logs_dir, exist_ok=True)

                shell(
                    '''
                    prokka %s \
                    --force \
                    --centre X \
                    --compliant \
                    --cpus {threads} \
                    --outdir %s \
                    --locustag %s \
                    --prefix %s \
                    --kingdom {params.kingdom} \
                    2> %s
                    ''' % (bin_fa,
                           output_dir,
                           bin_id, bin_id,
                           log_file))

                if os.path.exists(gff_file):
                    gff_count += 1

            if gff_count == len(bin_list):
                shell('''touch {output.done}''')


    rule predict_bins_gene_prokka_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene/{{assembler}}.{{binner}}.prokka.out/{sample}/done"),
                sample=SAMPLES.index.unique())
        output:
            html = os.path.join(
                config["output"]["predict"],
                "report/bins_gene_{assembler}.{binner}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "report/bins_gene_{assembler}.{binner}.multiqc.out/prokka_multiqc_report_data"))
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/report/bins_gene_{assembler}.{binner}.multiqc.prokka.log")
        params:
            input_dir = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prokka.out/"),
            output_dir = os.path.join(
                config["output"]["predict"],
                "report/bins_gene_{assembler}.{binner}.multiqc.out")
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
                    "bins_gene/{assembler}.{binner}.prokka.out/{sample}/done"),
                os.path.join(
                    config["output"]["predict"],
                    "report/bins_gene_{assembler}.{binner}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "report/bins_gene_{assembler}.{binner}.multiqc.out/prokka_multiqc_report_data")],
                   assembler=ASSEMBLERS,
                   binner=BINNERS,
                   sample=SAMPLES.index.unique())

else:
    rule predict_bins_gene_prokka_all:
        input:


rule predict_scaftigs_gene_all:
    input:
        rules.predict_scaftigs_gene_prodigal_all.input,
        rules.predict_scaftigs_gene_prokka_all.input,

        rules.assembly_all.input


rule predict_bins_gene_all:
    input:
        rules.predict_bins_gene_prodigal_all.input,
        rules.predict_bins_gene_prokka_all.input,

        rules.binning_all.input


rule predict_all:
    input:
        rules.predict_scaftigs_gene_all.input,
        rules.predict_bins_gene_all.input
