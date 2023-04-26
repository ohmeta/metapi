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
        params:
            phamb_utils = config["params"]["annotation"]["dbscan_swa"]["phamb_utils"],
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
            -m {params.min_binsize}
            '''


    rule annotation_prophage_dbscan_swa:
        input:
            os.path.join(
                config["output"]["annotation"],
                "dbscan_swa/{binning_group}.{assembler}.mags/vamb_bins.{batchid}.fna")
        output:
            done = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.output/vamb_bins.{batchid}/done"),
            fna = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/{batchid}.fna"),
            faa = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/{batchid}.faa"),
            summary = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/{batchid}.summary.txt")
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

            cat {params.outdir}/*.fna > {output.fna}
            cat {params.outdir}/*.faa > {output.faa}

            head -1 {params.outdir}/*1*summary.txt > {output.summary}
            tail -n +2 -q {params.outdir}/*prophage_summary.txt >> {output.summary}

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
            os.path.join(
                config["output"]["annotation"],
                "dbscan_swa/{binning_group}.{assembler}.output/done")
        shell:
            '''
            touch {output}
            '''


    rule annotation_prophage_dbscan_swa_all:
        input:
            expand(os.path.join(
                config["output"]["annotation"],
                "dbscan_swa/{binning_group}.{assembler}.output/done"),
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS)

else:
    rule annotation_prophage_dbscan_swa_all:
        input:


rule annotation_all:
    input:
        rules.annotation_prophage_dbscan_swa_all.input


localrules:
    annotation_prophage_dbscan_swa_all,
    annotation_all
