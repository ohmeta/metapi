if config["params"]["predict"]["prodigal"]["do"]:
    rule predict_scaftigs_gene:
        input:
            os.path.join(config["output"]["assembly"],
                         "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz")
        output:
            pep = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.pep.faa"),
            cds = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.cds.ffn"),
            gff = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.cds.gff"),
            score = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.score.gff")
        log:
            os.path.join(config["output"]["predict"],
                         "logs/scaftigs_gene/{sample}.{assembler}.prodigal.log")
        params:
            format = config["params"]["predict"]["prodigal"]["format"],
            mode = config["params"]["predict"]["prodigal"]["mode"]
        shell:
            '''
            zcat {input} | \
            prodigal \
            -a {output.pep} \
            -d {output.cds} \
            -o {output.gff} \
            -s {output.score} \
            -f {params.format} \
            -p {params.mode} -q \
            2> {log}
            '''


    rule predict_scaftigs_gene_all:
        input:
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.pep.faa"),
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.cds.ffn"),
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.cds.gff"),
                os.path.join(
                    config["output"]["predict"],
                    "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.score.gff")],
                   assembler=ASSEMBLERS,
                   sample=SAMPLES.index.unique())

else:
    rule predict_scaftigs_gene_all:
        input:


if config["params"]["predict"]["prokka"]["do"]:
    rule predict_bins_gene:
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
                                    "logs/bins_gene/{sample}"),
            kingdom = config["params"]["predict"]["prokka"]["kingdom"],
            metagenome = "--metagenome" \
                if config["params"]["predict"]["prokka"]["metagenome"] \
                   else ""
        threads:
            config["params"]["predict"]["prokka"]["threads"]
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
                    {params.metagenome} \
                    2> %s
                    ''' % (bin_fa,
                           output_dir,
                           bin_id, bin_id,
                           log_file))

                if os.path.exists(gff_file):
                    gff_count += 1

            if gff_count == len(bin_list):
                shell('''touch {output.done}''')


    rule predict_bins_gene_multiqc:
        input:
            expand(
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene/{{assembler}}.{{binner}}.prokka.out/{sample}/done"),
                sample=SAMPLES.index.unique())
        output:
            html = os.path.join(
                config["output"]["predict"],
                "bins_gene_report/{assembler}.{binner}.multiqc.out/prokka_multiqc_report.html"),
            data_dir = directory(os.path.join(
                config["output"]["predict"],
                "bins_gene_report/{assembler}.{binner}.multiqc.out/prokka_multiqc_report_data"))
        log:
            os.path.join(
                config["output"]["predict"],
                "logs/bins_gene_report/{assembler}.{binner}.multiqc.prokka.log")
        params:
            input_dir = os.path.join(
                config["output"]["predict"],
                "bins_gene/{assembler}.{binner}.prokka.out/"),
            output_dir = os.path.join(
                config["output"]["predict"],
                "bins_gene_report/{assembler}.{binner}.multiqc.out")
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


    rule predict_bins_gene_all:
        input:
            expand([
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene/{assembler}.{binner}.prokka.out/{sample}/done"),
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene_report/{assembler}.{binner}.multiqc.out/prokka_multiqc_report.html"),
                os.path.join(
                    config["output"]["predict"],
                    "bins_gene_report/{assembler}.{binner}.multiqc.out/prokka_multiqc_report_data")],
                   assembler=ASSEMBLERS,
                   binner=BINNERS,
                   sample=SAMPLES.index.unique())

else:
    rule predict_bins_gene_all:
        input:


rule predict_all:
    input:
        rules.predict_scaftigs_gene_all.input,
        rules.predict_bins_gene_all.input
