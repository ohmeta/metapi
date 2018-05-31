def _get_fastq(wildcards, units, read_pair="fq1"):
    return units.loc[(wildcards.sample, wildcards.unit), [read_pair]].dropna()[0]

rule simulate_reads:
    input:
        ref = config["simulate"]["ref"]
    output:
        r1 = lambda wildcards: _get_fastq(wildcards, units, "fq1"),
        r2 = lambda wildcards: _get_fastq(wildcards, units, "fq2")
    shell:
        '''
        randomreads.sh ref={input.ref} len=100 illuminanames addslash paired reads=10000000 metagenome
        '''
