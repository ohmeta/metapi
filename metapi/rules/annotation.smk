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


    checkpoint annotation_prophage_dbscan_swa_merge:
        input:
            aggregate_dbscan_swa_output
        output:
            fna = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/prophage.fna"),
            faa = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/prophage.faa"),
            summary = os.path.join(config["output"]["annotation"], "dbscan_swa/{binning_group}.{assembler}.prophage/prophage_summary.tsv")
        shell:
            '''
            for outdone in {input}
            do
                outdir=$(dirname $outdone)
                summary=$outdir/test_DBSCAN-SWA_prophage_summary.txt
                head -1 $summary > {output.summary}
                break
            done

            for outdone in {input}
            do
                outdir=$(dirname $outdone)

                FNA=$outdir/test_DBSCAN-SWA_prophage.fna
                FAA=$outdir/test_DBSCAN-SWA_prophage.faa
                summary=$outdir/test_DBSCAN-SWA_prophage_summary.txt

                [ -s $FNA ] && cat $FNA >> {output.fna}
                [ -s $FAA ] && cat $FAA >> {output.faa}
                [ -s $summary ] && tail -n +2 -q $summary >> {output.summary}
            done
            '''

    def get_dbscan_swa_merged_output(wildcards):
         checkpoint_output = checkpoints.annotation_prophage_dbscan_swa_merge.get(**wildcards).output.fna
         return expand(os.path.join(
             config["output"]["annotation"], 
             "dbscan_swa/{binning_group}.{assembler}.prophage/prophage.fna"),
             binning_group=wildcards.binning_group,
             assembler=wildcards.assembler)

    checkpoint annotation_prophage_dbscan_swa_distribute:
        input:
            metadata = os.path.join(config["output"]["assembly"],
                 "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz"),
            all_fna = get_dbscan_swa_merged_output
        output:
            os.path.join(config["output"]["identify"], "vmags/{binning_group}.{assembly_group}.{assembler}/dbscan_swa/distribution_done")
        params:
            working_dir = os.path.join(config["output"]["identify"], "vmags/{binning_group}.{assembly_group}.{assembler}/dbscan_swa"),
            assembly_fna = os.path.join(config["output"]["identify"], "vmags/{binning_group}.{assembly_group}.{assembler}/dbscan_swa/{binning_group}.{assembly_group}.{assembler}.dbscan_swa.combined.fa"),
            assembly_group = "{assembly_group}"
        run:
            shell("rm -rf {params.working_dir}")
            shell("mkdir -p {params.working_dir}")
            shell("touch {params.assembly_fna}")

            import pandas as pd
            from Bio import SeqIO

            ### record assembly_group : alias ###
            tab = pd.read_table(input.metadata)
            assembly_vamb_id = {vamb_id : binning_assembly.split(".")[-1] for vamb_id, binning_assembly in zip(tab.iloc[:,1], tab.iloc[:, 0])}

            ### read the prophage fna ###
            for record in SeqIO.parse(input.all_fna[0], 'fasta'):
                desc = record.description
                vamb_id = desc.split("|")[0].split("C")[0]
                if assembly_vamb_id[vamb_id] != params.assembly_group:
                    # print(vamb_id, assembly_vamb_id[vamb_id])
                    continue
                with open(params.assembly_fna, "a") as f:
                    f.write(record.format("fasta"))

            shell("gzip -f {params.assembly_fna}")
            shell("touch {output}")

    rule annotation_prophage_dbscan_swa_all:
        input:
            expand(os.path.join(
                config["output"]["identify"],
                "vmags/{binning_group}.{assembly_group}.{assembler}/dbscan_swa/distribution_done"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST,
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
