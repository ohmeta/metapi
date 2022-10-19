if config["params"]["checkv"]["do"]:
    rule checkv_download_db:
        output:
            expand([
                os.path.join(config["params"]["checkv"]["db"], "genome_db/checkv_{gfile}"),
                os.path.join(config["params"]["checkv"]["db"], "hmm_db/{hfile}")],
                gfile=["error.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.tsv", "source.tsv"],
                hfile=["checkv_hmms.tsv", "genome_lengths.tsv"])
        params:
            db = config["params"]["checkv"]["db"]
        benchmark:
            os.path.join(config["output"]["check"], "benchmark/checkv/checkv_download_db.benchmark.txt")
        log:
            os.path.join(config["output"]["check"], "logs/checkv/checkv_download_db.log")
        conda:
            config["envs"]["checkv"]
        shell:
            '''
            checkv download_database $(dirname {params.db}) \
            > {log} 2>&1
            '''
    

    rule checkv:
        input:
            db = expand([
                os.path.join(config["params"]["checkv"]["db"], "genome_db/checkv_{gfile}"),
                os.path.join(config["params"]["checkv"]["db"], "hmm_db/{hfile}")],
                gfile=["error.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.tsv", "source.tsv"],
                hfile=["checkv_hmms.tsv", "genome_lengths.tsv"]),
            viral = os.path.join(
                config["output"]["identify"],
                "vmags/{binning_group}.{assembly_group}.{assembler}/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.combined.fa")
        output:
            os.path.join(config["output"]["check"],
                         "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/checkv_done")
        benchmark:
            os.path.join(config["output"]["check"],
                         "benchmark/checkv/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.checkv.benchmark.txt")
        log:
            os.path.join(config["output"]["check"],
                         "logs/checkv/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.checkv.log")
        params:
            db = config["params"]["checkv"]["db"],
            outdir = os.path.join(config["output"]["check"], "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}")
        conda:
            config["envs"]["checkv"]
        threads:
            config["params"]["checkv"]["threads"]
        shell:
            '''
            rm -rf {params.outdir}

            checkv end_to_end \
            {input.viral} \
            {params.outdir} \
            -t {threads} \
            -d {params.db} \
            >{log} 2>&1

            touch {output}
            '''


    rule checkv_postprocess:
        input:
            done = os.path.join(config["output"]["check"],
                                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/checkv_done")
        output:
            vmag = os.path.join(config["output"]["check"],
                                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/vMAG_hmq.fa")
        params:
            binning_group = "{binning_group}",
            assembly_group = "{assembly_group}",
            assembler = "{assembler}",
            identifier = "{identifier}"
        run:
            import os
            import subprocess
            import pandas as pd
            from Bio import SeqIO

            label = f'''{params.binning_group}.{params.assembly_group}.{params.assembler}.{params.identifier}'''
            contig_num = 0

            checkv_dir = os.path.dirname(input.done)
            quality_summary_f = os.path.join(checkv_dir, "quality_summary.tsv")
            proviruses_f = os.path.join(checkv_dir, "proviruses.fna")
            viruses_f = os.path.join(checkv_dir, "viruses.fna")

            if os.path.exists(quality_summary_f):
                quality_summary = pd.read_csv(quality_summary_f, sep="\t")
                hmq_contig_ids = list(quality_summary.query('checkv_quality=="High-quality" or checkv_quality=="Medium-quality"')["contig_id"])
                print(f'''Identified {len(hmq_contig_ids)} high or medium quality vMAGs''')

                quality_summary = quality_summary.set_index("contig_id")

                if len(hmq_contig_ids) > 0:
                    with open(output.vmag, "w") as oh:
                        if os.path.exists(proviruses_f):
                            print(f'''Found {proviruses_f}''')
                            rc_proviruses = SeqIO.index_db(":memory:", proviruses_f, "fasta")
                            print(f'''Loaded {proviruses_f}''')
                            for contig_id in hmq_contig_ids:
                                if contig_id in rc_proviruses:
                                    contig_num += 1
                                    quality = quality_summary.loc[contig_id, "checkv_quality"]
                                    rc = rc_proviruses[contig_id]
                                    rc.id = f'''{label}.{contig_num}|proviruses|{quality}'''
                                    SeqIO.write(rc, oh, "fasta")
                        if os.path.exists(viruses_f):
                            print(f'''Found {viruses_f}''')
                            rc_viruses = SeqIO.index_db(":memory:", viruses_f, "fasta")
                            print(f'''Loaded {viruses_f}''')
                            for contig_id in hmq_contig_ids:
                                if contig_id in rc_viruses:
                                    contig_num += 1
                                    quality = quality_summary.loc[contig_id, "checkv_quality"]
                                    rc = rc_viruses[contig_id]
                                    rc.id = f'''{label}.{contig_num}|viruses|{quality}'''
                                    SeqIO.write(rc, oh, "fasta")
                    print(contig_num)
                else:
                    subprocess.run(f'''touch {output.vmag}''', shell=True)
            else:
                subprocess.run(f'''touch {output.vmag}''', shell=True)


    checkv_df_list = []
    for identifier in config["params"]["checkv"]["checkv_identifier"]:
        checkv_df = ASSEMBLY_GROUPS.copy()
        checkv_df["identifier"] = identifier 
        checkv_df_list.append(checkv_df)
    CHECKV_GROUPS = pd.concat(checkv_df_list, axis=0)


    rule checkv_all:
        input:
            expand([
                os.path.join(config["output"]["check"],
                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/checkv_done"),
                os.path.join(config["output"]["check"],
                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/vMAG_hmq.fa")],
                zip,
                binning_group=CHECKV_GROUPS["binning_group"],
                assembly_group=CHECKV_GROUPS["assembly_group"],
                assembler=CHECKV_GROUPS["assembler"],
                identifier=CHECKV_GROUPS["identifier"])
 
else:
    rule checkv_all:
        input:


localrules:
    checkv_download_db