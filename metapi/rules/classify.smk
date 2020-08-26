if config["params"]["classify"]["kraken2"]["do"]:
    rule classify_short_reads_kraken2:
        input:
            reads = assembly_input_with_short_reads,
            database = expand(os.path.join(
                config["params"]["classify"]["kraken2"]["database"], "{db}"),
                              db = ["hash.k2d", "taxo.k2d", "opts.k2d"])
        output:
            table = protected(os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.table.gz")),
            report = protected(os.path.join(
                config["output"]["classify"],
                "short_reads/{sample}.kraken2.out/{sample}.kraken2.report.gz"))
        log:
            os.path.join(config["output"]["classify"],
                         "logs/{sample}.kraken2.log")
        params:
            paired = "--paired" if IS_PE else "",
            database = config["params"]["classify"]["kraken2"]["database"],
            quick = "--quick" \
                if config["params"]["classify"]["kraken2"]["quick"] \
                   else "",
            memory_mapping = "--memory-mapping" \
                if config["params"]["classify"]["kraken2"]["memory_mapping"] \
                   else "",
            use_names = "--use-names" \
                if config["params"]["classify"]["kraken2"]["use_names"] \
                   else "",
            use_mpa_style = "--use-mpa-style" \
                if config["params"]["classify"]["kraken2"]["use_mpa_style"] \
                   else "",
            report_zero_counts = "--report-zero-counts" \
                if config["params"]["classify"]["kraken2"]["report_zero_counts"] \
                   else "",
            confidence = config["params"]["classify"]["kraken2"]["confidence"],
            min_base_quality = config["params"]["classify"]["kraken2"]["min_base_quality"],
            min_hit_groups = config["params"]["classify"]["kraken2"]["min_hit_groups"],
            unclassified_out = "--unclassified-out %s" % \
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.unclassified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["classify"]["kraken2"]["unclassified_out"] \
                       else "",
            classified_out = "--classified-out %s" % \
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.classified%s.fq" \
                    % "#" if IS_PE else "") \
                    if config["params"]["classify"]["kraken2"]["classified_out"] \
                       else ""
        threads:
            config["params"]["classify"]["threads"]
        run:
            import os

            shell(
                '''
                kraken2 \
                {params.quick} \
                {params.memory_mapping} \
                {params.use_mpa_style} \
                {params.use_names} \
                {params.report_zero_counts} \
                --threads {threads} \
                --db {params.database} \
                --confidence {params.confidence} \
                --minimum-base-quality {params.min_base_quality} \
                --minimum-hit-groups {params.min_hit_groups} \
                {params.unclassified_out} \
                {params.classified_out} \
                --output %s \
                --report %s \
                --gzip-compressed \
                {params.paired} \
                {input.reads} \
                2> {log}
                ''' % (os.path.splitext(output.table)[0],
                       os.path.splitext(output.report)[0]))

            output_dir = os.path.dirname(output.report)
            with os.scandir(output_dir) as itr:
                for entry in itr:
                    if entry.is_file():
                        shell('''pigz -p {threads} %s''' \
                              % os.path.join(output_dir, entry.name))


    rule classify_short_reads_kraken2_all:
        input:
            expand([
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.table.gz"),
                os.path.join(
                    config["output"]["classify"],
                    "short_reads/{sample}.kraken2.out/{sample}.kraken2.report.gz")],
                   sample=SAMPLES.index.unique()),

            rules.rmhost_all.input,
            rules.qcreport_all.input

else:
    rule classify_short_reads_kraken2_all:
        input:


if config["params"]["classify"]["gtdbtk"]["do"]:
    checkpoint classify_hmq_bins_gtdbtk_prepare:
        input:
            bins_hmq = os.path.join(
                config["output"]["checkm"],
                "report/{assembler}_{binner_checkm}_bins_hmq.tsv")
        output:
            out_dir = directory(
                os.path.join(config["output"]["classify"],
                             "bins_hmq_input/{assembler}_{binner_checkm}"))
        params:
            batch_num = config["params"]["classify"]["gtdbtk"]["batch_num"]
        run:
            import pandas as pd
            import os

            df = pd.read_csv(input.bins_hmq, names=["path"])
            df["id"] = df.apply(lambda x: os.path.basename(x["path"]), axis=1)

            os.makedirs(output.out_dir, exist_ok=True)

            if len(df) > 0:
                for batch_id in range(0, len(df), params.batch_num):
                    df_ = df.iloc[batch_id:batch_id + params.batch_num][["path", "id"]]
                    df_.to_csv(os.path.join(output.out_dir, "bins_hmq_%d.tsv" % batch_id),
                               sep='\t', index=False, header=None)
            else:
                shell('''touch {output.out_dir}/bins_hmq_0.tsv''')


    rule classify_hmq_bins_gtdbtk:
        input:
            bins_hmq = os.path.join(
                config["output"]["classify"],
                "bins_hmq_input/{assembler}_{binner_checkm}/bins_hmq_{batchid}.tsv")
        output:
            done = os.path.join(
                config["output"]["classify"],
                "table/{assembler}.{binner_checkm}.gtdbtk.out.{batchid}/done")
        wildcard_constraints:
            batchid="\d+"
        log:
            os.path.join(
                config["output"]["classify"],
                "logs/bins_hmq_{batchid}.{assembler}.{binner_checkm}.gtdbtk.log")
        params:
            bin_suffix = "fa",
            out_dir = os.path.join(
                config["output"]["classify"],
                "table/{assembler}.{binner_checkm}.gtdbtk.out.{batchid}")
        threads:
            config["params"]["classify"]["threads"]
        shell:
            '''
            gtdbtk classify_wf \
            --batchfile {input.bins_hmq}/ \
            --out_dir {params.out_dir} \
            --extension {params.bin_suffix} \
            --cpus {threads} \
            > {log}

            touch {output.done}
            '''


    def aggregate_gtdbtk_report_input(wildcards):
        checkpoint_output = checkpoints.classify_hmq_bins_gtdbtk_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["classify"],
            "table/{assembler}.{binner_checkm}.gtdbtk.out.{batchid}/done"),
                      assembler=wildcards.assembler,
                      binner_checkm=wildcards.binner_checkm,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(
                                                os.path.join(checkpoint_output,
                                                             "bins_hmq_{batchid}.tsv")).batchid])))


    rule classify_hmq_bins_gtdbtk_report:
        input:
            aggregate_gtdbtk_report_input
        output:
            table_gtdb_a = os.path.join(
                config["output"]["classify"],
                "report/bins_hmq.{assembler}.{binner_checkm}.gtdbtk.archaea.gtdb.tsv"),
            table_gtdb_b = os.path.join(
                config["output"]["classify"],
                "report/bins_hmq.{assembler}.{binner_checkm}.gtdbtk.bacteria.gtdb.tsv")
        threads:
            8
        run:
            import os

            ar122_list = []
            bac120_list = []

            for i in input:
                ar122_tsv = os.path.join(os.path.dirname(i), "gtdbtk.ar122.summary.tsv")
                bac120_tsv = os.path.join(os.path.dirname(i), "gtdbtk.bac120.summary.tsv")
                if os.path.exists(ar122_tsv):
                    ar122_list.append(ar122_tsv)
                if os.path.exists(bac120_tsv):
                    bac120_list.append(bac120_tsv)

            metapi.merge(ar122_list, metapi.parse, threads, output=output.table_gtdb_a)
            metapi.merge(bac120_list, metapi.parse, threads, output=output.table_gtdb_b)


    rule single_classify_hmq_bins_gtdbtk_all:
        input:
            expand(
                os.path.join(
                    config["output"]["classify"],
                    "report/bins_hmq.{assembler}.{binner_checkm}.gtdbtk.{taxonomy}.{system}.tsv"),
                taxonomy=["archaea", "bacteria"],
                system=["gtdb"],
                assembler=ASSEMBLERS,
                binner_checkm=BINNERS_CHECKM),

            rules.checkm_all.input,

else:
    rule single_classify_hmq_bins_gtdbtk_all:
        input:


rule single_classify_all:
    input:
        rules.classify_short_reads_kraken2_all.input,
        rules.single_classify_hmq_bins_gtdbtk_all.input
