if config["params"]["taxonomic"]["gtdbtk"]["do"]:
    checkpoint taxonomic_gtdbtk_prepare:
        input:
            rep_genomes_info = os.path.join(config["output"]["dereplicate"],
                         "report/bacteriome/checkm_table_genomes_info.{assembler}.derep.tsv.gz")
        output:
            mags_dir = directory(os.path.join(config["output"]["taxonomic"], "mags_input/{assembler}"))
        params:
            batch_num = config["params"]["taxonomic"]["gtdbtk"]["batch_num"]
        run:
            metapi.gtdbtk_prepare(input.rep_genomes_info,
                                  params.batch_num,
                                  output.mags_dir)


    rule taxonomic_gtdbtk:
        input:
            mags_input = os.path.join(config["output"]["taxonomic"], "mags_input/{assembler}/mags_input.{batchid}.tsv"),
            gtdb_data_path = expand(os.path.join(
                config["params"]["taxonomic"]["gtdbtk"]["gtdb_data_path"], "{gtdbtk_dir}"),
                gtdbtk_dir = ["fastani", "markers", "masks", "metadata",
                              "mrca_red", "msa", "pplacer", "radii", "taxonomy"])
        output:
            os.path.join(config["output"]["taxonomic"], "table/gtdbtk/gtdbtk.out.{assembler}.{batchid}/gtdbtk_done")
        wildcard_constraints:
            batchid="\d+"
        conda:
            config["envs"]["gtdbtk"]
        log:
            os.path.join(config["output"]["taxonomic"], "logs/gtdbtk/gtdbtk.{assembler}.{batchid}.log")
        benchmark:
            os.path.join(config["output"]["taxonomic"], "benchmark/gtdbtk/gtdbtk.{assembler}.{batchid}.benchmark.txt")
        params:
            gtdb_data_path = config["params"]["taxonomic"]["gtdbtk"]["gtdb_data_path"],
            pplacer_threads = config["params"]["taxonomic"]["gtdbtk"]["pplacer_threads"]
        threads:
            config["params"]["taxonomic"]["threads"]
        shell:
            '''
            export GTDB_DATA_PATH={params.gtdb_data_path}

            outdir=$(dirname {output})

            gtdbtk classify_wf \
            --batchfile {input.mags_input} \
            --out_dir $outdir \
            --extension fa \
            --cpus {threads} \
            --pplacer_cpus {params.pplacer_threads} \
            > {log} 2>&1

            if [ -f $outdir/classify/gtdbtk.bac120.summary.tsv ];
            then
                pushd $outdir
                ln -s classify/gtdbtk.bac120.summary.tsv gtdbtk.bacteria.summary.tsv
                popd
            else
                pushd $outdir
                touch gtdbtk.bacteria.summary.tsv
                popd
            fi

            if [ -f $outdir/classify/gtdbtk.ar53.summary.tsv ];
            then
                pushd $outdir
                ln -s classify/gtdbtk.ar53.summary.tsv gtdbtk.archaea.summary.tsv
                popd
            elif [ -f $outdir/classify/gtdbtk.ar122.summary.tsv ];
            then
                pushd $outdir 
                ln -s classify/gtdbtk.ar122.summary.tsv gtdbtk.archaea.summary.tsv
                popd
            else
                pushd $outdir
                touch gtdbtk.archaea.summary.tsv
                popd
            fi

            touch {output}
            '''


    def aggregate_gtdbtk_report_input(wildcards):
        checkpoint_output = checkpoints.taxonomic_gtdbtk_prepare.get(**wildcards).output[0]

        return expand(os.path.join(
            config["output"]["taxonomic"],
            "table/gtdbtk/gtdbtk.out.{assembler}.{batchid}/gtdbtk_done"),
                      assembler=wildcards.assembler,
                      batchid=list(set([i.split("/")[0] \
                                        for i in glob_wildcards(
                                                os.path.join(checkpoint_output,
                                                             "mags_input.{batchid}.tsv")).batchid])))


    rule taxonomic_gtdbtk_report:
        input:
            gtdb_done = aggregate_gtdbtk_report_input,
            rep_genomes_info = os.path.join(config["output"]["dereplicate"],
                                            "report/bacteriome/checkm_table_genomes_info.{assembler}.derep.tsv.gz")
        output:
            table_gtdb = os.path.join(config["output"]["taxonomic"],
                                      "report/gtdbtk/MAGs_hmq.rep.{assembler}.gtdbtk.gtdb.tsv"),
            table_ncbi = os.path.join(config["output"]["taxonomic"],
                                      "report/gtdbtk/MAGs_hmq.rep.{assembler}.gtdbtk.ncbi.tsv"),
            table_all = os.path.join(config["output"]["taxonomic"],
                                     "report/gtdbtk/MAGs_hmq.rep.{assembler}.gtdbtk.all.tsv")
        params:
            metadata_archaea = config["params"]["taxonomic"]["gtdbtk"]["metadata_archaea"],
            metadata_bacteria = config["params"]["taxonomic"]["gtdbtk"]["metadata_bacteria"],
            gtdb_to_ncbi_script = config["params"]["taxonomic"]["gtdbtk"]["gtdb_to_ncbi_script"]
        conda:
            config["envs"]["gtdbtk"]
        threads:
            8
        script:
            "../wrappers/gtdbtk_postprocess.py"
 

    rule taxonomic_gtdbtk_all:
        input:
            expand(
                os.path.join(
                    config["output"]["taxonomic"],
                    "report/gtdbtk/MAGs_hmq.rep.{assembler}.gtdbtk.{system}.tsv"),
                system=["gtdb", "ncbi", "all"],
                assembler=ASSEMBLERS,
                binner_checkm=BINNERS_CHECKM),

            #rules.checkm_all.input,


    localrules:
        taxonomic_gtdbtk_prepare,
        taxonomic_gtdbtk_report

else:
    rule taxonomic_gtdbtk_all:
        input:


rule taxonomic_all:
    input:
        rules.taxonomic_gtdbtk_all.input


localrules:
    taxonomic_gtdbtk_all,
    taxonomic_all