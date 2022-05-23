MULTIALIGN_GROUP = SAMPLES.reset_index().loc[:, ["sample_id", "assembly_group", "binning_group"]].drop_duplicates()

MULTIBINING_INDEX = {}
for binning_group in SAMPLES_BINNING_GROUP_LIST:
    MULTIBINING_INDEX[binning_group] = {}
    assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, binning_group))
    count = 0
    for assembly_group in assembly_groups:
        count += 1 
        MULTIBINING_INDEX[binning_group][assembly_group] = f'''S{count}'''

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
            config["envs"]["minimap2"]
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
            config["envs"]["minimap2"]
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
            config["envs"]["minimap2"]
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
            config["envs"]["metabat2"]
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
            metapi.combine_jgi(input.jgi, output.matrix)


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
            cuda_module = config["params"]["binning"]["vamb"]["cuda_module"],
            use_cuda_module = int(config["params"]["binning"]["vamb"]["use_cuda_module"]),
            allow_small_scaftigs = 1 if config["params"]["binning"]["vamb"]["allow_small_scaftigs"] else 0,
            external_params = config["params"]["binning"]["vamb"]["external_params"]
        shell:
            '''
            set +e

            rm -rf {params.outdir}
            mkdir -p {params.outdir_base}

            nums=`zcat {input.scaftigs} | grep -c "^>"`

            if [ $nums -lt 4096 ];
            then
                echo "The total number of contigs of {input.scaftigs} is $nums, less than 4096" > {log} 2>&1
                echo "See here for help: https://github.com/RasmussenLab/vamb/issues/35" >> {log} 2>&1

                if [ {params.allow_small_scaftigs} -eq 0 ];
                then
                    mkdir -p {params.outdir}
                    touch {output.binning_done}
                    echo "Allow small scaftigs: False" >> {log} 2>&1
                    echo "Touch binning_done" >> {log} 2>&1
                    exit 0
                else
                    echo "Allow small scaftigs: True" >> {log} 2>&1
                    echo "Maybe you need to adjust the number of epochs and start batch size" >> {log} 2>&1
                    echo "Running vamb" >> {log} 2>&1
                fi
            else
                echo "The total number of contigs of {input.scaftigs} is $nums, greater than 4096" > {log} 2>&1
                echo "Running vamb" >> {log} 2>&1
            fi


            if [ {params.use_cuda_module} -eq 1 ];
            then
                module load {params.cuda_module}
                echo "module load {params.cuda_module}" >> {log} 2>&1
                which nvcc >> {log} 2>&1
            fi

            if [ "{params.cuda}" == "--cuda" ];
            then
                lspci | grep -i nvidia >> {log} 2>&1
                which python >> {log} 2>&1
                which vamb >> {log} 2>&1

                python -c 'import torch;print(torch.__file__)' >> {log} 2>&1
                python -c 'import torch;print(f"Torch CUDA: {{torch.cuda.is_available()}}")' >> {log} 2>&1
                python -c 'from torch.utils.cpp_extension import CUDA_HOME;print(CUDA_HOME)' >> {log} 2>&1
                python -c 'import os; print(os.environ.get("CUDA_PATH"))' >> {log} 2>&1
            fi

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
            >> {log} 2>&1

            echo "Running vamb completed" >> {log} 2>&1
            echo "Touch binning_done" >> {log} 2>&1
            touch {output.binning_done}

            exit 0
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
            import sys
            import pandas as pd

            metadata = []
            assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))
            assembly_index = int(assembly_groups.index(params.assembly_group)) + 1
            assembly_index = f'''S{assembly_index}'''

            ## Double check
            if assembly_index != MULTIBINING_INDEX[params.binning_group][params.assembly_group]:
                sys.exit("assembly_group index error")

            outdir = os.path.dirname(output.binning_done)
            bins_dir = os.path.dirname(input.binning_done)
            os.makedirs(outdir, exist_ok=True)
            bin_index = 0

            if os.path.exists(f'{bins_dir}/bins'):
                fna_list = sorted(glob(f'{bins_dir}/bins/{assembly_index}C*.fna'))

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