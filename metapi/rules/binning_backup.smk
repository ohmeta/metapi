"""
    checkpoint binning_vamb_postprocess_run:
        input:
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/{binning_group}.{assembler}.vamb.out/binning_done")
        output:
            success = os.path.join(config["output"]["multisplit_binning"],
                "bins/{binning_group}.{assembler}.vamb.out/postprocess_run"),
            metadata = os.path.join(
                config["output"]["multisplit_binning"],
                "bins/{binning_group}.{assembler}.vamb.out/cluster.metadata.tsv")
        benchmark:
            os.path.join(config["output"]["multisplit_binning"],
                         "benchmark/binning_vamb_postprocess.{binning_group}.{assembler}.benchmark.txt")
        params:
            bins_from = config["output"]["multisplit_binning"],
            bins_to = config["output"]["binning"],
            binning_group = "{binning_group}",
            assembler = "{assembler}"
        run:
            from glob import glob
            import os
            import pandas as pd

            metadata = []
            assembly_index = 0
            assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))

            for assembly_group in assembly_groups:
                outdir = os.path.join(params.bins_to, f'''bins/{assembly_group}.{params.assembler}.out/vamb''')
                os.makedirs(outdir, exist_ok=True)

                assembly_index += 1
                bin_index = 0
                bins_dir = os.path.dirname(input)
                fna_list = sorted(glob(f'{bins_dir}/bins/S{assembly_index}C*.fna'))

                # link method
                #for fna in fna_list:
                #    count_ += 1
                #    fna_source = f"../../../../../{fna}"
                #    fna_dist = f"{i}.{params.assembler}.vamb.bin.{count_}.fa"
                #    shell(
                #        f'''
                #        pushd {outdir} && \
                #        ln -s {fna_source} {fna_dist} && \
                #        popd
                #        ''')

                # copy method, rename 
                for fna in fna_list:
                    bin_index += 1
                    # bin_id = os.path.basename(fna).split(".")[0]
                    # bin_id = os.path.basename(fna).split(".")[0].split("C")[-1]
                    fna_dist = os.path.join(outdir, f"{assembly_group}.{params.assembler}.vamb.bin.{bin_index}.fa")
                    metadata.append((os.path.abspath(fna), os.path.abspath(fna_dist)))
                    shell(f'''cat {fna} | seqkit replace -p "^S\d+C" > {fna_dist}''')

                shell(f'''touch {outdir}/binning_done''')

            pd.DataFrame(metadata, columns=["vamb_bin", "vamb_postprocess_bin"])\
              .to_csv(output.metadata, sep='\t', index=False)

            shell('''touch {output.success}''')


    def aggregate_vamb_postprocess_output(wildcards):
        checkpoint_output = checkpoints.binning_vamb_postprocess_run.get(**wildcards).output.success
        binning_group = os.path.dirname(checkpoint_output).split(".")[0]

        return expand(os.path.join(
            config["output"]["binning"],
            "bins/{assembly_group}.{assembler}.out/vamb/binning_done"),
            assembler=wildcards.assembler,
            assembly_group=metapi.get_assembly_group_by_binning_group(SAMPLES, binning_group))

    
    rule binning_vamb_postprocess_done:
        input:
            aggregate_vamb_postprocess_output
        output:
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/{binning_group}.{assembler}.vamb.out/postprocess_done")
        shell:
            '''
            touch {output} 
            '''


    rule binning_vamb_all:
        input:
            expand([
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "bins/{binning_group}.{assembler}.vamb.out/{results}")],
                #os.path.join(
                #    config["output"]["binning"],
                #    "bins/{assembly_group}.{assembler}.out/vamb")],
                #assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST,
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                results=["clusters.tsv",
                         "latent.npz",
                         "lengths.npz",
                         "log.txt",
                         "model.pt",
                         "mask.npz",
                         "tnf.npz",
                         "bins",
                         "cluster.metadata.tsv",
                         "postprocess_success_run",
                         "postprocess_success_done"])

"""

"""
MULTIALIGN_GROUP = SAMPLES.reset_index().loc[:, ["sample_id", "assembly_group", "binning_group"]].drop_duplicates()

MULTIBINNING_GROUP = SAMPLES.reset_index().loc[:, ["assembly_group", "binning_group"]].drop_duplicates()

multibinning_df_list = []
for assembler in ASSEMBLERS:
    multibinning_df = MULTIBINNING_GROUP.copy()
    multibinning_df["assembler"] = assembler
    multibinning_df_list.append(multibinning_df)
MULTIBINNING_GROUPS = pd.concat(multibinning_df_list, axis=0)


    rule binning_vamb_postprocess:
        input:
            binning_done = os.path.join(config["output"]["multisplit_binning"],
                         "bins/{binning_group}.{assembler}.vamb.out/binning_done")
        output:
            metadata = os.path.join(config["output"]["multisplit_binning"],
                "bins/{binning_group}.{assembler}.vamb.out/bins_{assembly_group}/cluster.metadata.tsv"),
            binning_done = os.path.join(config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/vamb/binning_done")
        benchmark:
            os.path.join(config["output"]["multisplit_binning"],
                         "benchmark/binning_vamb_postprocess.{binning_group}.{assembler}/{assembly_group}.benchmark.txt")
        params:
            bins_from = config["output"]["multisplit_binning"],
            bins_to = config["output"]["binning"],
            binning_group = "{binning_group}",
            assembler = "{assembler}",
            assembly_group = "{assembly_group}"
        run:
            from glob import glob
            import os
            import pandas as pd

            metadata = []
            assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))
            assembly_index = assembly_groups.index(params.assembly_group)

            outdir = os.path.dirname(output.binning_done)
            bins_dir = os.path.dirname(input.binning_done)
            os.makedirs(outdir, exist_ok=True)
            bin_index = 0
            fna_list = sorted(glob(f'{bins_dir}/bins/S{assembly_index}C*.fna'))

            for fna in fna_list:
                bin_index += 1
                # bin_id = os.path.basename(fna).split(".")[0]
                # bin_id = os.path.basename(fna).split(".")[0].split("C")[-1]
                fna_dist = os.path.join(outdir, f'''{params.assembly_group}.{params.assembler}.vamb.bin.{bin_index}.fa''')
                metadata.append((os.path.abspath(fna), os.path.abspath(fna_dist)))
                shell(f'''cat {fna} | seqkit replace -p "^S\d+C" > {fna_dist}''')

            shell(f'''touch {output.binning_done}''')

            pd.DataFrame(metadata, columns=["vamb_bin", "vamb_postprocess_bin"])\
              .to_csv(output.metadata, sep='\t', index=False)


    rule binning_vamb_all:
        input:
            expand(
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "bins/{binning_group}.{assembler}.vamb.out/{results}"),
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                results=["clusters.tsv",
                         "latent.npz",
                         "lengths.npz",
                         "log.txt",
                         "model.pt",
                         "mask.npz",
                         "tnf.npz"]),
            expand([
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "bins/{binning_group}.{assembler}.vamb.out/bins_{assembly_group}/cluster.metadata.tsv"),
                os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}.{assembler}.out/vamb/binning_done")],
                zip,
                binning_group=MULTIBINNING_GROUPS["binning_group"],
                assembler=MULTIBINNING_GROUPS["assembler"],
                assembly_group=MULTIBINNING_GROUPS["assembly_group"])
"""