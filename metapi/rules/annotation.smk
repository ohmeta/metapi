if config["params"]["annotation"]["dbscan_swa"]["do"]:

    checkpoint annotation_prophage_dbscan_swa_prepare:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
            vamb_done = os.path.join(
                config["output"]["binning"],
                "mags_vamb/{binning_group}.{assembler}/binning_done")
        output:
            mags_dir = directory(os.path.join(config["output"]["annotation"],
                "dbscan_swa/{binning_group}.{assembler}.mags"))
        log:
            os.path.join(
                config["output"]["annotation"],
                "logs/dbscan_swa_prepare/{binning_group}.{assembler}.dbscan_swa_prepare.log")
        benchmark:
            os.path.join(
                config["output"]["annotation"],
                "benchmark/dbscan_swa_prepare/{binning_group}.{assembler}.dbscan_swa_prepare.benchmark.txt")
        params:
            phamb_utils = config["params"]["annotation"]["dbscan_swa"]["phamb_utils"],
            batch_num = config["params"]["annotation"]["dbscan_swa"]["batch_num"],
            min_binsize = config["params"]["annotation"]["dbscan_swa"]["min_binsize"],
            cluster_tsv = os.path.join(
                config["output"]["binning"],
                "mags_vamb/{binning_group}.{assembler}/clusters.tsv")
        shell:
            '''
            python {params.phamb_utils} \
            {input.scaftigs} \
            {params.cluster_tsv} \
            {output} \
            -m {params.min_binsize} \
            -b {params.batch_num} \
            >{log} 2>&1
            '''


    rule annotation_prophage_dbscan_swa:
        input:
            os.path.join(
                config["output"]["annotation"],
                "dbscan_swa/{binning_group}.{assembler}.mags/vamb_bins.{batchid}.fna")
        output:
            done = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.output/vamb_bins.{batchid}/done")
        log:
            os.path.join(config["output"]["annotation"], "logs/dbscan_swa/{binning_group}.{assembler}.{batchid}.dbscan_swa.log")
        benchmark:
            os.path.join(config["output"]["annotation"], "benchmark/dbscan_swa/{binning_group}.{assembler}.{batchid}.dbscan_swa.benchmark.txt")
        params:
            dbscan_swa_script = config["params"]["annotation"]["dbscan_swa"]["script"],
            outdir = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.output/vamb_bins.{batchid}"),
            prefix = "test"
        threads:
            config["params"]["annotation"]["threads"]
        shell:
            '''
            python {params.dbscan_swa_script} \
            --input {input} \
            --output {params.outdir} \
            --thread_num {threads} \
            --prefix {params.prefix} \
            --add_annotation none \
            >{log} 2>&1

            touch {output.done}
            '''


    def aggregate_dbscan_swa_output(wildcards):
        checkpoint_output = checkpoints.annotation_prophage_dbscan_swa_prepare.get(**wildcards).output.mags_dir

        return expand(os.path.join(
            config["output"]["annotation"],
            "dbscan_swa/{binning_group}.{assembler}.output/vamb_bins.{batchid}/done"),
            binning_group=wildcards.binning_group,
            assembler=wildcards.assembler,
            batchid=list(set([i for i in glob_wildcards(os.path.join(checkpoint_output, "vamb_bins.{batchid}.fna")).batchid])))


    rule annotation_prophage_dbscan_swa_merge:
        input:
            aggregate_dbscan_swa_output
        output:
            fna = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/prophage.fna"),
            faa = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/prophage.faa"),
            summary = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/prophage_summary.tsv")
        shell:
            '''
            outdone={input}[0]
            outdir=$(dirname $outdone)
            summary=$outdir/test_DBSCAN-SWA_prophage_summary.txt
            head -1 $summary > {output.summary}

            for outdone in {input}
            do
                outdir=$(dirname $outdone)

                FNA=$outdir/test_DBSCAN-SWA_prophage.fna
                FAA=$outdir/test_DBSCAN-SWA_prophage.faa
                summary=$outdir/test_DBSCAN-SWA_prophage_summary.txt

                [ -s $fna ] && cat $FNA >> {output.fna}
                [ -s $faa ] && cat $FAA >> {output.faa}
                [ -s $summary ] && tail -n +2 -q $summary >> {output.summary}
            done
            '''


    rule annotation_prophage_dbscan_swa_all:
        input:
            expand(os.path.join(
                config["output"]["annotation"],
                "dbscan_swa/{binning_group}.{assembler}.prophage/{results}"),
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                results=["prophage.fna", "prophage.faa", "prophage_summary.tsv"])

else:
    rule annotation_prophage_dbscan_swa_all:
        input:


rule annotation_all:
    input:
        rules.annotation_prophage_dbscan_swa_all.input


localrules:
    annotation_prophage_dbscan_swa_all,
    annotation_all
