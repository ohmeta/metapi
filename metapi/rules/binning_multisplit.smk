MULTIALIGN_GROUP = SAMPLES.reset_index().loc[:, ["sample_id", "assembly_group", "binning_group"]].drop_duplicates()

MULTIBINNING_GROUP = SAMPLES.reset_index().loc[:, ["assembly_group", "binning_group"]].drop_duplicates()

multibinning_df_list = []
for assembler in ASSEMBLERS:
    multibinning_df = MULTIBINNING_GROUP.copy()
    multibinning_df["assembler"] = assembler
    multibinning_df_list.append(multibinning_df)
MULTIBINNING_GROUPS = pd.concat(multibinning_df_list, axis=0)


def get_reads_for_multisplit_binning(wildcards, step):
    samples_id_list = metapi.get_samples_id_by_binning_group(SAMPLES, wildcards.binning_group)
    short_reads = get_short_reads_list(step, samples_id_list)
    return short_reads


def multisplit_binning_input_with_short_reads(wildcards):
    if RMHOST_DO:
        return get_reads_for_multisplit_binning(wildcards, "rmhost", False, False)
    elif TRIMMING_DO:
        return get_reads_for_multisplit_binning(wildcards, "trimming", False, False)
    else:
        return get_reads_for_multisplit_binning(wildcards, "raw", False, False)


# reference: https://github.com/RasmussenLab/vamb/blob/master/workflow/vamb.snake.conda.py
if config["params"]["binning"]["vamb"]["do"]:
    rule binning_vamb_combine_scaftigs:
        input:
            lambda wildcards: expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{{assembler}}.out/{assembly_group}.{{assembler}}.scaftigs.fa.gz"),
                   assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            os.path.join(
                config["output"]["multisplit_binning"],
                "scaftigs/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.fa.gz")
        benchmark:
            os.path.join(
                config["output"]["multisplit_binning"],
                "benchmark/binning_vamb_combine_scaftigs.{binning_group}.{assembler}.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(
                config["output"]["multisplit_binning"],
                "logs/binning_vamb_combine_scaftigs.{binning_group}.{assembler}.log")
        params:
            min_contig = config["params"]["binning"]["vamb"]["min_contig"]
        shell:
            '''
            concatenate.py {output} {input} -m {params.min_contig} 2> {log}
            '''


    rule binning_vamb_dict_scaftigs:
        input:
            os.path.join(
                config["output"]["multisplit_binning"],
                "scaftigs/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.fa.gz")
        output:
            os.path.join(
                config["output"]["multisplit_binning"],
                "index/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.dict")
        benchmark:
            os.path.join(
                config["output"]["multisplit_binning"],
                "benchmark/binning_vamb_dict_scaftigs.{binning_group}.{assembler}.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/binning_vamb_dict_scaftigs.{binning_group}.{assembler}.log")
        shell:
            '''
            samtools dict {input} | cut -f1-3 > {output} 2> {log}
            '''


    rule binning_vamb_index_scaftigs:
        input:
            os.path.join(
                config["output"]["multisplit_binning"],
                "scaftigs/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.fa.gz")
        output:
            os.path.join(
                config["output"]["multisplit_binning"],
                "index/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.minimap2.mmi")
        benchmark:
            os.path.join(
                config["output"]["multisplit_binning"],
                "benchmark/binning_vamb_index_scaftigs.{binning_group}.{assembler}.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/binning_vamb_index_scaftigs.{binning_group}.{assembler}.log")
        params:
            index_size = config["params"]["binning"]["vamb"]["index_size"]
        shell:
            '''
            minimap2 -I {params.index_size} -d {output} {input} 2> {log}
            '''


    rule binning_vamb_align_scaftigs:
        input:
            reads = alignment_input_with_short_reads,
            scaftigs_index = os.path.join(
                config["output"]["multisplit_binning"],
                "index/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.minimap2.mmi"),
            scaftigs_dict = os.path.join(
                config["output"]["multisplit_binning"],
                "index/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.dict")
        output:
            flagstat = os.path.join(
                config["output"]["multisplit_binning"],
                "report/flagstat_minimap2/{binning_group}.{assembler}/{sample}.align2combined_scaftigs.flagstat"),
            bam = temp(os.path.join(config["output"]["multisplit_binning"],
                                    "bam/{binning_group}.{assembler}.out/{sample}/{sample}.align2combined_scaftigs.bam"))
        priority:
            28
        benchmark:
            os.path.join(config["output"]["multisplit_binning"],
                         "benchmark/minimap2/{binning_group}.{assembler}/{sample}.align2combined_scaftigs.benchmark.txt")
        params:
            bam_dir = os.path.join(config["output"]["multisplit_binning"], "bam/{binning_group}.{assembler}.out/{sample}") 
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/alignment/{binning_group}.{assembler}/{sample}.align2combined_scaftigs.log")
        threads:
            config["params"]["alignment"]["threads"]
        shell:
            '''
            rm -rf {params.bam_dir}
            mkdir -p {params.bam_dir}

            minimap2 -t {threads} -ax sr {input.scaftigs_index} {input.reads} -N 5 2> {log} |
            tee >(samtools flagstat \
                  -@{threads} - > {output.flagstat}) | \
            grep -v "^@" | \
            cat {input.scaftigs_dict} - | \
            samtools view -F 3584 -b - > {output.bam} 2> {log}
            
            #samtools sort -m 3G -@{threads} -T {output.bam} -O BAM -o {output.bam} -
            '''


    rule binning_vamb_sort_bam:
        input:
            bam = os.path.join(
                config["output"]["multisplit_binning"],
                "bam/{binning_group}.{assembler}.out/{sample}/{sample}.align2combined_scaftigs.bam")
        output:
             bam = os.path.join(
                config["output"]["multisplit_binning"],
                "bam/{binning_group}.{assembler}.out/{sample}/{sample}.align2combined_scaftigs.sorted.bam") \
                if config["params"]["binning"]["vamb"]["save_bam"] else \
                   temp(os.path.join(
                       config["output"]["multisplit_binning"],
                       "bam/{binning_group}.{assembler}.out/{sample}/{sample}.align2combined_scaftigs.sorted.bam"))
        priority:
            29
        benchmark:
            os.path.join(config["output"]["multisplit_binning"],
                         "benchmark/minimap2/{binning_group}.{assembler}/{sample}.sort_bam.benchmark.txt")
        params:
            bam_dir = os.path.join(config["output"]["multisplit_binning"], "bam/{binning_group}.{assembler}.out/{sample}") 
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/alignment/{binning_group}.{assembler}/{sample}.sort_bam.log")
        threads:
            config["params"]["alignment"]["threads"]
 
        shell:
            '''
            rm -rf {params.bam_dir}/temp*
            
            samtools sort {input.bam} -m 3G -@{threads} -T {params.bam_dir}/temp -O BAM -o {output.bam} 2> {log}
            
            rm -rf {params.bam_dir}/temp*
            '''
 

    rule binning_vamb_align_scaftigs_report:
        input:
            expand(
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "report/flagstat_minimap2/{binning_group}.{{assembler}}/{sample}.align2combined_scaftigs.flagstat"),
                zip,
                binning_group=MULTIALIGN_GROUP["binning_group"],
                sample=MULTIALIGN_GROUP["sample_id"])
        output:
            flagstat = os.path.join(config["output"]["multisplit_binning"],
                                    "report/alignment_flagstat_{assembler}.tsv")
        run:
            input_list = [str(i) for i in input]
            output_str = str(output)
            metapi.flagstats_summary(input_list, 2, output=output.flagstat)


    rule binning_vamb_coverage:
        input:
            bam = os.path.join(
                config["output"]["multisplit_binning"],
                "bam/{binning_group}.{assembler}.out/{sample}/{sample}.align2combined_scaftigs.sorted.bam")
        output:
            jgi = os.path.join(
                config["output"]["multisplit_binning"],
                "coverage/{binning_group}.{assembler}.out/{sample}.align2combined_scaftigs.jgi")
        priority:
            30
        benchmark:
            os.path.join(
                config["output"]["multisplit_binning"],
                "benchmark/jgi_summarize_bam_contig_depths/{binning_group}.{assembler}/{sample}.jgi_summarize_bam_contig_depths.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/coverage/{binning_group}.{assembler}/{sample}.align2combined_scaftigs.jgi.coverage.log")
        shell:
            '''
            jgi_summarize_bam_contig_depths \
            --noIntraDepthVariance --outputDepth {output.jgi} {input.bam} 2> {log}
            '''


    rule binning_vamb_gen_abundance_matrix:
        input:
            jgi = lambda wildcards: expand(os.path.join(
                config["output"]["multisplit_binning"],
                "coverage/{{binning_group}}.{{assembler}}.out/{sample}.align2combined_scaftigs.jgi"),
                sample=metapi.get_samples_id_by_binning_group(SAMPLES, wildcards.binning_group))
        output:
            matrix = os.path.join(config["output"]["multisplit_binning"],
                                  "matrix/{binning_group}.{assembler}.align2combined_scaftigs.jgi.abundance.matrix.tsv")
        benchmark:
            os.path.join(
                config["output"]["multisplit_binning"],
                "benchmark/binning_vamb_gen_abundance_matrix.{binning_group}.{assembler}.benchmark.txt")
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/binning_vamb_gen_abundance_matrix.{binning_group}.{assembler}.log")
        run:
            import pandas as pd
            import numpy as np

            first = False
            #jgi_df_list = []
            #for jgi in input.jgi:
            #    if not first:
            #        # jgi format
            #        # contigName\tcontigLen\ttotalAvgDepth\t{sample_id}.align2combined_scaftigs.sorted.bam
            #        jgi_df_first = pd.read_csv(jgi, sep="\t")\
            #                     .loc[:, ["contigName", "contigLen", "totalAvgDepth"]]\
            #                     .dtype({"contigName": str, "contigLen": np.int32, "totalAvgDepth": np.float32})\
            #                     .set_index("contigName")
            #        jgi_df = pd.read_csv(jgi, sep="\t").iloc[:, [0, 3]]\
            #                   .dtype({"contigName": str})
            #        jgi_df[jgi_df.columns[1]] = jgi_df[jgi_df.columns[1]].astype(np.float32)
            #        jgi_df_list = [jgi_df_first, jgi_df.set_index("contigName")]
            #        first = True
            #    else:
            #        jgi_df = pd.read_csv(jgi, sep="\t").iloc[:, [0, 3]]\
            #                   .dtype({"contigName": str})
            #        jgi_df[jgi_df.columns[1]] = jgi_df[jgi_df.columns[1]].astype(np.float32)
            #        jgi_df_list.append(jgi_df.set_index("contigName"))
            ## big table, huge memory
            #pd.concat(jgi_df_list, axis=1).reset_index().to_csv(output.matrix, sep="\t", index=False)

            matrix_list = []
            for jgi in input.jgi:
                if not first:
                    first = True
                    with open(jgi, 'r') as ih:
                        for line in ih:
                            line_list = line.strip().split("\t")
                            matrix_list.append(line_list)
                else:
                    with open(jgi, 'r') as ih:
                        count = -1
                        for line in ih:
                            count += 1
                            line_list = line.strip().split("\t")
                            matrix_list[count].append(line_list[3])

            with open(output.matrix, 'w') as oh:
                for i in matrix_list:
                    oh.write("\t".join(i) + "\n")


    rule binning_vamb_prepare_all:
        input:
            expand([
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "scaftigs/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "matrix/{binning_group}.{assembler}.align2combined_scaftigs.jgi.abundance.matrix.tsv"),
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "report/alignment_flagstat_{assembler}.tsv")],
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS)


    rule binning_vamb:
        input:
            scaftigs = os.path.join(config["output"]["multisplit_binning"],
                "scaftigs/{binning_group}.{assembler}.combined.out/{binning_group}.{assembler}.combined.scaftigs.fa.gz"),
            matrix = os.path.join(config["output"]["multisplit_binning"],
                "matrix/{binning_group}.{assembler}.align2combined_scaftigs.jgi.abundance.matrix.tsv")
        output:
            #vamb_output = expand(
            #    os.path.join(config["output"]["multisplit_binning"],
            #                 "bins/{{binning_group}}.{{assembler}}.vamb.out/{results}"),
            #    results=["clusters.tsv", "latent.npz", "lengths.npz",
            #             "log.txt", "model.pt", "mask.npz", "tnf.npz"]),
            binning_done = os.path.join(config["output"]["multisplit_binning"],
                                        "bins/{binning_group}.{assembler}.vamb.out/binning_done")
        benchmark:
            os.path.join(config["output"]["multisplit_binning"],
                         "benchmark/vamb/{binning_group}.{assembler}.vamb.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/binning/{binning_group}.{assembler}.vamb.binning.log")
        threads:
            config["params"]["binning"]["threads"]
        params:
            outdir = os.path.join(config["output"]["multisplit_binning"],
                "bins/{binning_group}.{assembler}.vamb.out"),
            outdir_base = os.path.join(config["output"]["multisplit_binning"], "bins"),
            min_contig = config["params"]["binning"]["vamb"]["min_contig"],
            min_fasta = config["params"]["binning"]["vamb"]["min_fasta"],
            cuda = "--cuda" if config["params"]["binning"]["vamb"]["cuda"] else "",
            external_params = config["params"]["binning"]["vamb"]["external_params"]
        shell:
            '''
            rm -rf {params.outdir}
            mkdir -p {params.outdir_base}

            if [[ `zcat {input.scaftigs} | grep -c "^>"` -lt 4096 ]];
            then
                mkdir -p {params.outdir}
                touch {output.binning_done}
            else
                vamb \
                {params.cuda} \
                -p {threads} \
                --outdir {params.outdir} \
                --fasta {input.scaftigs} \
                --jgi {input.matrix} \
                -o C \
                -m {params.min_contig} \
                --minfasta {params.min_fasta} \
                {params.external_params} \
                2> {log}
                touch {output.binning_done}
            fi
            '''


    rule binning_vamb_postprocess:
        input:
            binning_done = os.path.join(config["output"]["multisplit_binning"],
                         "bins/{binning_group}.{assembler}.vamb.out/binning_done")
        output:
            metadata = os.path.join(config["output"]["multisplit_binning"],
                "bins/{binning_group}.{assembler}.vamb.out/bins_{assembly_group}/cluster.metadata.tsv"),
            binning_done = os.path.join(config["output"]["binning"],
                "bins/{assembly_group}__{binning_group}.{assembler}.out/vamb/binning_done")
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
            assembly_index = int(assembly_groups.index(params.assembly_group)) + 1

            outdir = os.path.dirname(output.binning_done)
            bins_dir = os.path.dirname(input.binning_done)
            os.makedirs(outdir, exist_ok=True)
            bin_index = 0

            if os.path.exists(f'{bins_dir}/bins'):
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
                results=[
                    #"clusters.tsv",
                    #"latent.npz",
                    #"lengths.npz",
                    #"log.txt",
                    #"model.pt",
                    #"mask.npz",
                    #"tnf.npz",
                    "binning_done"]),
            expand([
                os.path.join(
                    config["output"]["multisplit_binning"],
                    "bins/{binning_group}.{assembler}.vamb.out/bins_{assembly_group}/cluster.metadata.tsv"),
                os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}__{binning_group}.{assembler}.out/vamb/binning_done")],
                zip,
                binning_group=MULTIBINNING_GROUPS["binning_group"],
                assembler=MULTIBINNING_GROUPS["assembler"],
                assembly_group=MULTIBINNING_GROUPS["assembly_group"])

else:
    rule binning_vamb_prepare_all:
        input:

    rule binning_vamb_all:
        input:


localrules:
    binning_vamb_prepare_all,
    binning_vamb_postprocess,
    binning_vamb_all