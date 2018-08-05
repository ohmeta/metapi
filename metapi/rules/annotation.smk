rule prokka_bins:
    input:
        directory(os.path.join(config["results"]["binning"]["bins"], "{sample}.metabat2_out"))
    output:
        directory(os.path.join(config["results"]["annotation"]["prokka"], "{sample}.prokka_out"))
    params:
        logdir = config["logs"]["annotation"]["prokka"],
        kingdom = config["params"]["annotation"]["prokka"]["kingdom"],
        metagenome = "--metagenome" if config["params"]["annotation"]["prokka"]["metagenome"] else ""
    threads:
        config["params"]["annotation"]["prokka"]["threads"]
    run:
        import glob
        import os
        import sys
        import subprocess
        import shutil
        bin_list = glob.glob(input[0] + "/*bin*fa")
        for bin in bin_list:
            bin_id = os.path.basename(bin.strip()).rstrip(".fa")
            log_file= os.path.join(params.logdir, bin_id + ".prokka.log")
            os.makedirs(params.logdir, exist_ok=True)
            cmd = "prokka %s --outdir %s --locustag %s --prefix %s --kingdom %s %s --cpus %d 2> %s" % (bin.strip(), output[0], bin_id, bin_id, params.kingdom, params.metagenome, threads, log_file)
            print(cmd)
            subprocess.Popen(cmd, shell=True)
