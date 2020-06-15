STEPS = ["raw"]
if TRIMMING_DO:
    STEPS += ["trimming"]
if RMHOST_DO:
    STEPS += ["rmhost"]

if config["params"]["qcreport"]["do"]:
    rule qcreport_summary:
        input:
            expand(os.path.join(config["output"]["qcreport"],
                                "{step}_stats.tsv"),
                   step=STEPS)
        output:
            os.path.join(config["output"]["qcreport"], "qc_stats.tsv")
        priority:
            30
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            df = metapi.merge(input, metapi.parse, threads)
            metapi.compute_host_rate(df, output=output[0])


    rule qcreport_plot:
        input:
            os.path.join(config["output"]["qcreport"], "qc_stats.tsv")
        output:
            os.path.join(config["output"]["qcreport"], "qc_reads_num_barplot.pdf")
        priority:
            30
        run:
            metapi.qc_bar_plot(input[0], "seaborn", output=output[0])


    rule qcreport_all:
        input:
            os.path.join(config["output"]["qcreport"], "qc_stats.tsv"),
            os.path.join(config["output"]["qcreport"], "qc_reads_num_barplot.tsv")

            rules.rmhost_all.input

else:
    rule qcreport_all:
        input:
