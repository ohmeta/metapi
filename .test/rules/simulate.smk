rule simulate_reads:
    input:
        ref = config["simulate"]["ref"]
    output:
        default = os.path.join(config["logs"]["simulate"], "done"),
        r1 = config["simulate"]["r1"],
        r2 = config["simulate"]["r2"]
    params:
        fa = os.path.basename(config["simulate"]["ref"]),
        dir = config["results"]["simulate"]
    shell:
        '''
        cd {params.dir}
        randomreads.sh ref={params.fa} len=100 illuminanames addslash paired reads=10000000 metagenome
        cd ../../
        touch {output.default}
        '''
