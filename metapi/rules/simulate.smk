if config["params"]["simulate"]["do"]:
    rule simulate_short_reads:
        input:
            genomes = lambda wildcards: metapi.get_simulate_info(SAMPLES, wildcards, "genome")
        output:
            r1 = os.path.join(config["output"]["simulate"],
                              "short_reads/{sample}.simulate.1.fq.gz"),
            r2 = os.path.join(config["output"]["simulate"],
                              "short_reads/{sample}.simulate.2.fq.gz"),
            abunf = os.path.join(config["output"]["simulate"],
                                 "abundance/{sample}.simulate.abundance.txt")
        log:
            os.path.join(config["output"]["simulate"], "logs/{sample}.iss.log")
        benchmark:
            os.path.join(config["output"]["simulate"], "benchmark/iss/{sample}.iss.benchmark.txt")
        conda:
            config["envs"]["simulate"]
        params:
            output_prefix = os.path.join(config["output"]["simulate"],
                                         "short_reads/{sample}"),
            model = lambda wildcards: metapi.get_simulate_info(SAMPLES, wildcards, "model")[0],
            reads_num = lambda wildcards: metapi.get_simulate_info(SAMPLES, wildcards, "reads_num")[0],
            abundance = lambda wildcards: metapi.get_simulate_info(SAMPLES, wildcards, "abundance")
        threads:
            config["params"]["simulate"]["threads"]
        script:
            "../wrappers/simulate_reads.py"


    rule simulate_all:
        input:
            expand([
                os.path.join(config["output"]["simulate"],
                             "short_reads/{sample}.simulate.{read}.fq.gz"),
                os.path.join(config["output"]["simulate"],
                             "abundance/{sample}.simulate.abundance.txt")],
                   read=["1", "2"],
                   sample=SAMPLES.index.unique())

else:
    rule simulate_all:
        input:


localrules:
    simulate_all