rule prokka_bins:
    input:
        os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.metabat2_out")
    output:
        default = os.path.join(config["results"]["annotation"]["prokka"], "{sample}.{assembler}.prokka_out/done")
    params:
        outdir = directory(os.path.join(config["results"]["annotation"]["prokka"], "{sample}.{assembler}.prokka_out")),
        logdir = config["logs"]["annotation"]["prokka"],
        kingdom = config["params"]["annotation"]["prokka"]["kingdom"],
        metagenome = "--metagenome" if config["params"]["annotation"]["prokka"]["metagenome"] else ""
    threads:
        config["params"]["annotation"]["prokka"]["threads"]
    run:
        import glob
        import os
        import time
        import subprocess
        bin_list = glob.glob(input[0] + "/*bin*fa")

        cmd_list = []
        for bin in bin_list:
            bin_id = os.path.basename(bin.strip()).rstrip(".fa")
            prokka_dir = os.path.join(params.outdir, bin_id)
            log_file= os.path.join(params.logdir, bin_id + ".prokka.log")
            os.makedirs(params.logdir, exist_ok=True)
            cmd = "prokka %s --outdir %s --locustag %s --prefix %s --kingdom %s %s --cpus %d 2> %s" % (bin.strip(), prokka_dir, bin_id, bin_id, params.kingdom, params.metagenome, threads, log_file)
            print(cmd)
            cmd_list.append(cmd)

        rc = subprocess.check_call("&&".join(cmd_list), shell=True)
        time.sleep(60)

        if rc == 0:
            with open(output.default, 'w') as out:
                out.write("Hello, Prokka!")


rule multiqc_prokka_bins:
    input:
        expand("{prokka}/{sample}.{assembler}.prokka_out/done",
               prokka=config["results"]["annotation"]["prokka"],
               assembler=config["params"]["assembler"],
               sample=_samples.index)
    output:
        html = os.path.join(config["results"]["annotation"]["multiqc_prokka"], "prokka_multiqc_report.html"),
        data_dir = directory(os.path.join(config["results"]["annotation"]["multiqc_prokka"], "prokka_multiqc_report_data"))
    log:
        os.path.join(config["logs"]["annotation"]["prokka"], "multiqc_prokka.log")
    params:
        inputdir = config["results"]["annotation"]["prokka"],
        outdir = config["results"]["annotation"]["multiqc_prokka"]
    shell:
        '''
        multiqc --cl_config "prokka_fn_snames: True" --outdir {params.outdir} --title prokka --module prokka {params.inputdir} 2> {log}
        '''
