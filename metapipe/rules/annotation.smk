_bins = parse_bins(config["results"]["binning"]["bins"])

rule prokka_bins:
    input:
        lambda wildcards: get_bin_id(_bins, wildcards, "path")
    output:
        file = expand("{prokka}/{{bin}}/{{bin}}.{suffix}",
                      prokka=config["results"]["annotation"]["prokka"],
                      suffix=["gff", "gbk", "fna", "faa", "ffn", "sqn", "fsa", "tbl", "err", "log", "txt", "tsv"])
    log:
        os.path.join(config["logs"]["annotation"]["prokka"], "{bin}.prokka.log")
    params:
        outdir = directory(os.path.join(config["results"]["annotation"]["prokka"], "{bin}")),
        kingdom = config["params"]["annotation"]["prokka"]["kingdom"],
        metagenome = "--metagenome" if config["params"]["annotation"]["prokka"]["metagenome"] else ""
    threads:
        config["params"]["annotation"]["prokka"]["threads"]
    shell:
        '''
        prokka {input} --outdir {params.outdir} --prefix {bin} --kingdom {params.kingdom} {params.metagenome} --cpus {threads} 2> {log}
        '''
