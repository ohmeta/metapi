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

        
    localrules:
        checkv_download_db


    rule checkv:
        input:
            db = expand([
                os.path.join(config["params"]["checkv"]["db"], "genome_db/checkv_{gfile}"),
                os.path.join(config["params"]["checkv"]["db"], "hmm_db/{hfile}")],
                gfile=["error.tsv", "reps.dmnd", "reps.faa", "reps.fna", "reps.tsv", "source.tsv"],
                hfile=["checkv_hmms.tsv", "genome_lengths.tsv"]),
            viral = os.path.join(
                config["output"]["identify"],
                "vmags/{binning_group}.{assembly_group}.{assembler}/{identifier}/{binning_group}.{assembly_group}.{assembler}.{identifier}.combined.fa.gz")
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

            FILESIZE=$(stat -c %s {input.viral})

            if [ $FILESIZE -gt 0 ];
            then
                checkv end_to_end \
                {input.viral} \
                {params.outdir} \
                -t {threads} \
                -d {params.db} \
                >{log} 2>&1
            fi

            touch {output}
            '''


    rule checkv_postprocess:
        input:
            done = os.path.join(config["output"]["check"],
                                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/checkv_done")
        output:
            vmag = os.path.join(config["output"]["check"],
                                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/vMAG_hmq.fa.gz")
        params:
            binning_group = "{binning_group}",
            assembly_group = "{assembly_group}",
            assembler = "{assembler}",
            identifier = "{identifier}"
        run:
            import os
            import gzip
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

                proviruses_df = quality_summary\
                .query('provirus=="Yes"')\
                .query('checkv_quality=="Complete" or checkv_quality=="High-quality" or checkv_quality=="Medium-quality"')\
                .set_index("contig_id")

                viruses_df = quality_summary\
                .query('provirus=="No"')\
                .query('checkv_quality=="Complete" or checkv_quality=="High-quality" or checkv_quality=="Medium-quality"')\
                .set_index("contig_id")

                print(f'''Identified {len(proviruses_df) + len(viruses_df)} complete, high or medium quality vMAGs''')

                with gzip.open(output.vmag, "wt") as oh:
                    if os.path.exists(proviruses_f):
                        proviruses_rc = SeqIO.index_db(":memory:", proviruses_f, "fasta")
                        for contig_id in proviruses_df.index.unique():
                            quality = proviruses_df.loc[contig_id, "checkv_quality"]
                            if contig_id in proviruses_rc:
                                contig_num += 1
                                rc = proviruses_rc[contig_id]
                                rc.id = f'''{label}.{contig_num}|proviruses|{quality}'''
                                SeqIO.write(rc, oh, "fasta")
                            else:
                                contig_id = f'''{contig_id}_1'''
                                if contig_id in proviruses_rc:
                                    contig_num += 1
                                    rc = proviruses_rc[contig_id]
                                    rc.id = f'''{label}.{contig_num}|proviruses|{quality}'''
                                    SeqIO.write(rc, oh, "fasta")
                                else:
                                    print(f'''Proviruses contig_id {contig_id} can't be found in {proviruses_f}, please check it!''')
                                    sys.exit(1)
                    else:
                        print(f'''No {proviruses_f} can be found''')

                    if os.path.exists(viruses_f):
                        viruses_rc = SeqIO.index_db(":memory:", viruses_f, "fasta")
                        for contig_id in viruses_df.index.unique():
                            quality = viruses_df.loc[contig_id, "checkv_quality"]
                            if contig_id in viruses_rc:
                                contig_num += 1
                                rc = viruses_rc[contig_id]
                                rc.id = f'''{label}.{contig_num}|viruses|{quality}'''
                                SeqIO.write(rc, oh, "fasta")
                            else:
                                print(f'''Viruses contig_id {contig_id} can't be found in {viruses_f}, please check it!''')
                                sys.exit(1)
                    else:
                        print(f'''No {viruses_f} can be found''')

            else:
                vmag = os.path.splitext(output.vmag)[0]
                subprocess.run(f'''touch {vmag}''', shell=True)
                subprocess.run(f'''pigz -f {vmag}''', shell=True)


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
                "data/checkv/{binning_group}.{assembly_group}.{assembler}/{identifier}/vMAG_hmq.fa.gz")],
                zip,
                binning_group=CHECKV_GROUPS["binning_group"],
                assembly_group=CHECKV_GROUPS["assembly_group"],
                assembler=CHECKV_GROUPS["assembler"],
                identifier=CHECKV_GROUPS["identifier"])
 
else:
    rule checkv_all:
        input:


rule check_all:
    input:
        rules.checkm_all.input,
        rules.checkv_all.input


localrules:
    check_all