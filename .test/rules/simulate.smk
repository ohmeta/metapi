rule simulate_reads:
    input:
        ref = config["simulate"]["ref"]
    output:
        default = os.path.join(config["logs"]["simulate"], "done"),
        r1 = config["simulate"]["r1"],
        r2 = config["simulate"]["r2"]
    shell:
        '''
        randomreads.sh ref={input.ref} len=100 illuminanames addslash paired reads=10000000 metagenome
        touch {output.default}
        '''
