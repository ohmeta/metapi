STEPS = ["raw"]
if TRIMMING_DO:
    STEPS += ["trimming"]
if RMHOST_DO:
    STEPS += ["rmhost"]

"""
if config["params"]["qcreport"]["do"]:
    rule qcreport:
        input:
            lambda wildcards: get_reads(wildcards, wildcards.step)
        output:
            os.path.join(config["output"]["{{step}}"],
                         "report/stats/{sample}_{step}_stats.tsv")
        params:
            sample_id = "{sample}",
            step = "{step}",
            fq_encoding = config["params"]["fq_encoding"]
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            shell(
                '''
                seqkit stats \
                --all \
                --basename \
                --tabular \
                --fq-encoding {params.fq_encoding} \
                --out-file {output} \
                --threads {threads} \
                {intput}
                ''')

            if IS_PE:
                metapi.change(output[0], params.sample_id, params.step,
                              "pe", ["fq1", "fq2"])
            else:
                metapi.change(output[0], params.sample_id, params.step,
                              "se", ["fq1"])


    rule qcreport_merge:
        input:
            expand(
                os.path.join(config["output"]["{{step}}"],
                             "report/stats/{sample}_{{step}}_stats.tsv"),
                sample=SAMPLES.index.unique())
        output:
            os.path.join(config["output"]["qcreport"], "{step}_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            metapi.merge(input, metapi.parse, threads,
                         save=True, output=output[0])


    rule qcreport_summary:
        input:
            expand(os.path.join(config["output"]["qcreport"],
                                "{step}_stats.tsv"),
                   step=STEPS)
        output:
            os.path.join(config["output"]["qcreport"], "qc_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            df = metapi.merge(input, metapi.parse, threads)
            metapi.compute_host_rate(df, save=True, output=output[0])


    rule qcreport_all:
        input:
            expand([
                os.path.join(config["output"]["qcreport"],
                             "{step}_stats.tsv"),
                os.path.join(config["output"]["qcreport"],
                             "qc_stats.tsv")],
                   step=STEPS)
"""
#else:
rule qcreport_all:
    input:
