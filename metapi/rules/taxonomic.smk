if config["params"]["taxonomic"]["gtdbtk"]["do"]:
    checkpoint taxonomic_gtdbtk_prepare:
        input:
            rep_genomes_info = os.path.join(config["output"]["dereplicate"],
                         "report/checkm_table_genomes_info.derep.tsv")
        output:
            bins_dir = directory(os.path.join(config["output"]["taxonomic"], "bins_input"))
        params:
            batch_num = config["params"]["taxonomic"]["gtdbtk"]["batch_num"]
        run:
            metapi.gtdbtk_prepare(input.rep_genomes_info, params.batch_num, output.bins_dir)


    rule taxonomic_gtdbtk:
        input:
            bins_input = os.path.join(config["output"]["taxonomic"], "bins_input/bins_input.{batchid}.tsv"),
            gtdb_data_path = expand(os.path.join(
                config["params"]["taxonomic"]["gtdbtk"]["gtdb_data_path"], "{gtdbtk_dir}"),
                gtdbtk_dir = ["fastani", "markers", "masks", "metadata",
                              "mrca_red", "msa", "pplacer", "radii", "taxonomy"])
        output:
            gtdbtk_done = os.path.join(config["output"]["taxonomic"], "table/gtdbtk.out.{batchid}/gtdbtk_done")
        wildcard_constraints:
            batchid="\d+"
        conda:
            config["envs"]["gtdbtk"]
        log:
            os.path.join(config["output"]["taxonomic"], "logs/gtdbtk.{batchid}.log")
        benchmark:
            os.path.join(config["output"]["taxonomic"], "benchmark/gtdbtk.{batchid}.benchmark.txt")
        params:
            out_dir = os.path.join(config["output"]["taxonomic"], "table/gtdbtk.out.{batchid}"),
            gtdb_data_path = config["params"]["taxonomic"]["gtdbtk"]["gtdb_data_path"],
            pplacer_threads = config["params"]["taxonomic"]["gtdbtk"]["pplacer_threads"]
        threads:
            config["params"]["taxonomic"]["threads"]
        shell:
            '''
            export GTDB_DATA_PATH={params.gtdb_data_path}

            gtdbtk classify_wf \
            --batchfile {input.bins_input} \
            --out_dir {params.out_dir} \
            --extension fa \
            --cpus {threads} \
            --pplacer_cpus {params.pplacer_threads} \
            2> {log} 2>&1

            touch {output.gtdbtk_done}
            '''


    def aggregate_gtdbtk_report_input(wildcards):
        checkpoint_output = checkpoints.taxonomic_gtdbtk_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["taxonomic"],
            "table/gtdbtk.out.{batchid}/gtdbtk_done"),
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(
                                                os.path.join(checkpoint_output,
                                                             "bins_input.{batchid}.tsv")).batchid])))


    rule taxonomic_gtdbtk_report:
        input:
            gtdb_table = aggregate_gtdbtk_report_input,
            rep_genomes_info = os.path.join(config["output"]["dereplicate"],
                                            "report/checkm_table_genomes_info.derep.tsv")
        output:
            table_gtdb = os.path.join(config["output"]["taxonomic"],
                                      "report/bins.hmq.rep.gtdbtk.gtdb.tsv"),
            table_ncbi = os.path.join(config["output"]["taxonomic"],
                                      "report/bins.hmq.rep.gtdbtk.ncbi.tsv"),
            table_all = os.path.join(config["output"]["taxonomic"],
                                     "report/bins.hmq.rep.gtdbtk.all.tsv")
        params:
            ar122_metadata = config["params"]["taxonomic"]["gtdbtk"]["ar122_metadata"],
            bac120_metadata = config["params"]["taxonomic"]["gtdbtk"]["bac120_metadata"],
            gtdb_to_ncbi_script = config["params"]["taxonomic"]["gtdbtk"]["gtdb_to_ncbi_script"]
        threads:
            8
        run:
            import os
            import pandas as pd
           
            gtdb_list = []
            ncbi_list = []

            for i in input.gtdb_table:
                out_dir = os.path.dirname(i)
                ar122_tsv = os.path.join(out_dir, "gtdbtk.ar122.summary.tsv")
                bac120_tsv = os.path.join(out_dir, "gtdbtk.bac120.summary.tsv")
           
                if os.path.exists(ar122_tsv):
                    gtdb_list.append(ar122_tsv)
                if os.path.exists(bac120_tsv):
                    gtdb_list.append(bac120_tsv)
           
                gtdb_to_ncbi_summary = os.path.join(out_dir, "gtdbtk.ncbi.summary.tsv")
                gtdb_to_ncbi_log = os.path.join(out_dir, "gtdbtk.to.ncbi.log")
           
                shell(
                    f'''
                    python {params.gtdb_to_ncbi_script} \
                    --gtdbtk_output_dir {out_dir} \
                    --output_file {gtdb_to_ncbi_summary} \
                    --ar122_metadata_file {params.ar122_metadata} \
                    --bac120_metadata_file {params.bac120_metadata} \
                    > {gtdb_to_ncbi_log}
                    ''')
           
                if os.path.exists(gtdb_to_ncbi_summary):
                    ncbi_list.append(gtdb_to_ncbi_summary)
           
            metapi.merge(gtdb_list, metapi.parse, threads, output=output.table_gtdb)
            metapi.merge(ncbi_list, metapi.parse, threads, output=output.table_ncbi)
           
            table_gtdb = pd.read_csv(output.table_gtdb, sep="\t").rename(columns={"classification": "GTDB classification"})
            table_ncbi = pd.read_csv(output.table_ncbi, sep="\t")
            table_rep_genomes_info = pd.read_csv(input.rep_genomes_info, sep="\t").rename(columns={"genome": "user_genome"})

            table_all = pd.merge(table_gtdb, table_ncbi, how="inner", on=["user_genome", "GTDB classification"]).\
                           merge(table_rep_genomes_info, how="inner", on="user_genome")
           
            table_all.to_csv(output.table_all, sep="\t", index=False)
           

    rule taxonomic_gtdbtk_all:
        input:
            expand(
                os.path.join(
                    config["output"]["taxonomic"],
                    "report/bins.hmq.rep.gtdbtk.{system}.tsv"),
                system=["gtdb", "ncbi", "all"],
                assembler=ASSEMBLERS,
                binner_checkm=BINNERS_CHECKM),

            #rules.checkm_all.input,

else:
    rule taxonomic_gtdbtk_all:
        input:


rule taxonomic_all:
    input:
        rules.taxonomic_gtdbtk_all.input


localrules:
    taxonomic_gtdbtk_prepare,
    taxonomic_gtdbtk_all,
    taxonomic_gtdbtk_report,
    taxonomic_all