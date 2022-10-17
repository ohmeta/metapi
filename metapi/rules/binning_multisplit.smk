MULTIBINING_INDEX = {}
for binning_group in SAMPLES_BINNING_GROUP_LIST:
    MULTIBINING_INDEX[binning_group] = {}
    assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, binning_group))
    count = 0
    for assembly_group in assembly_groups:
        count += 1 
        MULTIBINING_INDEX[binning_group][assembly_group] = f'''S{count}'''


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
                "scaftigs/{{binning_group}}.{assembly_group}.{{assembler}}/{{binning_group}}.{assembly_group}.{{assembler}}.scaftigs.fa.gz"),
                   assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/vamb/index/vamb_combine_scaftigs.{binning_group}.{assembler}.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(
                config["output"]["binning"],
                "logs/vamb/index/vamb_combine_scaftigs.{binning_group}.{assembler}.log")
        params:
            min_contig = config["params"]["binning"]["vamb"]["min_contig"]
        shell:
            '''
            concatenate.py {output} {input} -m {params.min_contig} 2> {log}
            '''


    rule binning_vamb_gen_metadata:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
        output:
            os.path.join(config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv")
        params:
            binning_group = "{binning_group}"
        run:
            import sys

            assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))

            with open(output[0], 'w') as oh:
                oh.write("binning_assembly_group\tvamb_id\n")
                count = 0
                for assembly_group in assembly_groups:
                    count += 1
                    vamb_id = f'''S{count}'''
                    # double check
                    if vamb_id != MULTIBINING_INDEX[params.binning_group][assembly_group]:
                        print("VAMB: sample id issue!")
                        sys.exit(1)
                    else:
                        oh.write(f'''{params.binning_group}.{assembly_group}\tS{count}\n''')


    rule binning_vamb_dict_scaftigs:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
        output:
            os.path.join(
                config["output"]["alignment"],
                "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.dict")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/vamb/index/samtools_dict_merged_scaftigs.{binning_group}.{assembler}.benchmark.txt")
        conda:
            config["envs"]["align"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/index/samtools_dict_merged_scaftigs.{binning_group}.{assembler}.log")
        shell:
            '''
            samtools dict {input} | cut -f1-3 > {output} 2> {log}
            '''


    rule binning_vamb_index_scaftigs:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
        output:
            os.path.join(
                config["output"]["alignment"],
                "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.minimap2.mmi")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/vamb/index/minimap2_index_merged_scaftigs.{binning_group}.{assembler}.benchmark.txt")
        conda:
            config["envs"]["align"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/index/minimap2_index_merged_scaftigs.{binning_group}.{assembler}.log")
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
                config["output"]["alignment"],
                "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.minimap2.mmi"),
            scaftigs_dict = os.path.join(
                config["output"]["alignment"],
                "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.dict")
        output:
            flagstat = os.path.join(
                config["output"]["alignment"],
                "report/flagstat_minimap2/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.flagstat"),
            bam = temp(os.path.join(config["output"]["alignment"],
                                    "bam_merged/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.bam"))
        priority:
            28
        params:
            bam_dir = os.path.join(config["output"]["alignment"],
                                   "bam_merged/{binning_group}.{assembler}"),
            sample = "{sample}"
        conda:
            config["envs"]["align"]
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/vamb/minimap2/{binning_group}.{assembler}.{sample}.align2merged_scaftigs.benchmark.txt")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/minimap2/{binning_group}.{assembler}.{sample}.align2merged_scaftigs.log")
        threads:
            config["params"]["alignment"]["threads"]
        shell:
            '''
            rm -rf {params.bam_dir}/{params.sample}.align2merged_scaftigs.bam*

            minimap2 -t {threads} -ax sr {input.scaftigs_index} {input.reads} -N 5 2> {log} |
            tee >(samtools flagstat \
                  -@{threads} - > {output.flagstat}) | \
            grep -v "^@" | \
            cat {input.scaftigs_dict} - | \
            samtools view -F 3584 -b - > {output.bam} 2> {log}
            '''


    rule binning_vamb_sort_bam:
        input:
            bam = os.path.join(
                config["output"]["alignment"],
                "bam_merged/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.bam")
        output:
             bam = os.path.join(
                config["output"]["alignment"],
                "bam_merged/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.sorted.bam") \
                if config["params"]["binning"]["vamb"]["save_bam"] else \
                   temp(os.path.join(
                       config["output"]["alignment"],
                       "bam_merged/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.sorted.bam"))
        priority:
            29
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/vamb/samtools/{binning_group}.{assembler}.{sample}.sort_bam.benchmark.txt")
        params:
            bam_dir = os.path.join(config["output"]["alignment"],
                                   "bam_merged/{binning_group}.{assembler}"),
            sample = "{sample}"
        conda:
            config["envs"]["align"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/samtools/{binning_group}.{assembler}.{sample}.sort_bam.log")
        threads:
            config["params"]["alignment"]["threads"]
        shell:
            '''
            rm -rf {params.bam_dir}/{params.sample}.align2merged_scaftigs.bam.temp*
            
            samtools sort {input.bam} -m 3G -@{threads} -T {params.bam_dir}/{params.sample}.align2merged_scaftigs.bam.temp -O BAM -o {output.bam} 2> {log}
            
            rm -rf {params.bam_dir}/{params.sample}.align2merged_scaftigs.bam.temp*
            '''
 

    rule binning_vamb_align_scaftigs_report:
        input:
            expand(
                os.path.join(
                    config["output"]["alignment"],
                    "report/flagstat_minimap2/{binning_group}.{{assembler}}/{sample}.align2merged_scaftigs.flagstat"),
                zip,
                binning_group=ALIGNMENT_GROUP["binning_group"],
                sample=ALIGNMENT_GROUP["sample_id"])
        output:
            flagstat = os.path.join(config["output"]["alignment"],
                                    "report/alignment_flagstat_{assembler}_minimap2.tsv")
        run:
            input_list = [str(i) for i in input]
            output_str = str(output)
            metapi.flagstats_summary(input_list, 2, output=output.flagstat)


    rule binning_vamb_coverage:
        input:
            bam = os.path.join(
                config["output"]["alignment"],
                "bam_merged/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.sorted.bam")
        output:
            jgi = os.path.join(
                config["output"]["binning"],
                "coverage/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.jgi")
        priority:
            30
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/vamb/jgi/{binning_group}.{assembler}.{sample}.jgi_summarize_bam_contig_depths.benchmark.txt")
        conda:
            config["envs"]["metabat2"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/jgi/{binning_group}.{assembler}.{sample}.align2merged_scaftigs.jgi.coverage.log")
        shell:
            '''
            jgi_summarize_bam_contig_depths \
            --noIntraDepthVariance --outputDepth {output.jgi} {input.bam} 2> {log}
            '''


    rule binning_vamb_gen_abundance_matrix:
        input:
            jgi = lambda wildcards: expand(os.path.join(
                config["output"]["binning"],
                "coverage/{{binning_group}}.{{assembler}}/{sample}.align2merged_scaftigs.jgi"),
                sample=sorted(metapi.get_samples_id_by_binning_group(SAMPLES, wildcards.binning_group)))
        output:
            matrix = os.path.join(config["output"]["binning"],
                                  "matrix/{binning_group}.{assembler}.align2merged_scaftigs.jgi.abundance.matrix.tsv")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/vamb/matrix/binning_vamb_gen_abundance_matrix.{binning_group}.{assembler}.benchmark.txt")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/matrix/binning_vamb_gen_abundance_matrix.{binning_group}.{assembler}.log")
        run:
            metapi.combine_jgi(input.jgi, output.matrix)


    rule binning_vamb_prepare_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["binning"],
                    "matrix/{binning_group}.{assembler}.align2merged_scaftigs.jgi.abundance.matrix.tsv"),
                os.path.join(
                    config["output"]["alignment"],
                    "report/alignment_flagstat_{assembler}_minimap2.tsv")],
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS)


    rule binning_vamb:
        input:
            scaftigs = os.path.join(config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
            matrix = os.path.join(config["output"]["binning"],
                "matrix/{binning_group}.{assembler}.align2merged_scaftigs.jgi.abundance.matrix.tsv")
        output:
            binning_done = os.path.join(config["output"]["binning"],
                                        "mags_vamb/{binning_group}.{assembler}/binning_done")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/vamb/vamb/{binning_group}.{assembler}.vamb.benchmark.txt")
        conda:
            config["envs"]["vamb"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/vamb/vamb/{binning_group}.{assembler}.vamb.binning.log")
        threads:
            config["params"]["binning"]["threads"]
        params:
            outdir = os.path.join(config["output"]["binning"], "mags_vamb/{binning_group}.{assembler}"),
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
            mkdir -p $(dirname {params.outdir})

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
            metadata = os.path.join(config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv"),
            binning_done = os.path.join(config["output"]["binning"],
                                        "mags_vamb/{binning_group}.{assembler}/binning_done")
        output:
            metadata = os.path.join(config["output"]["binning"],
                                    "mags_vamb/{binning_group}.{assembler}/bins_{assembly_group}/cluster.metadata.tsv"),
            binning_done = os.path.join(config["output"]["binning"],
                                        "mags/{binning_group}.{assembly_group}.{assembler}/vamb/binning_done")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/vamb/postprocess/binning_vamb_postprocess.{binning_group}.{assembly_group}.{assembler}.benchmark.txt")
        params:
            binning_group = "{binning_group}",
            assembly_group = "{assembly_group}",
            assembler = "{assembler}",
        run:
            from glob import glob
            import os
            import sys
            import pandas as pd

            binning_assembly_metadata = pd.read_csv(input.metadata, sep="\t").set_index("binning_assembly_group")
            assembly_index = binning_assembly_metadata.loc[f'''{params.binning_group}.{params.assembly_group}''', "vamb_id"]

            #assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))
            #assembly_index = int(assembly_groups.index(params.assembly_group)) + 1
            #assembly_index = f'''S{assembly_index}'''
            ## Double check
            #if assembly_index != MULTIBINING_INDEX[params.binning_group][params.assembly_group]:
            #    sys.exit("assembly_group index error")

            metadata = []

            outdir = os.path.dirname(output.binning_done)
            mags_dir = os.path.dirname(input.binning_done)
            os.makedirs(outdir, exist_ok=True)
            bin_index = 0

            if os.path.exists(f'{mags_dir}/bins'):
                fna_list = sorted(glob(f'{mags_dir}/bins/{assembly_index}C*.fna'))

                for fna in fna_list:
                    bin_index += 1
                    # bin_id = os.path.basename(fna).split(".")[0]
                    # bin_id = os.path.basename(fna).split(".")[0].split("C")[-1]
                    fna_dist = os.path.join(outdir, f'''{params.binning_group}.{params.assembly_group}.{params.assembler}.vamb.bin.{bin_index}.fa''')
                    metadata.append((os.path.abspath(fna), os.path.abspath(fna_dist)))
                    shell(f'''cat {fna} | seqkit replace -p "^S\d+C" > {fna_dist}''')

            shell(f'''touch {output.binning_done}''')

            pd.DataFrame(metadata, columns=["vamb_bin", "vamb_postprocess_bin"])\
              .to_csv(output.metadata, sep='\t', index=False)


    rule binning_vamb_all:
        input:
            rules.binning_vamb_prepare_all.input,
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "mags_vamb/{binning_group}.{assembler}/{results}"),
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
                    config["output"]["binning"],
                    "mags_vamb/{binning_group}.{assembler}/bins_{assembly_group}/cluster.metadata.tsv"),
                os.path.join(
                    config["output"]["binning"],
                    "mags/{binning_group}.{assembly_group}.{assembler}/vamb/binning_done")],
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"])

else:
    rule binning_vamb_prepare_all:
        input:

    rule binning_vamb_all:
        input:


localrules:
    binning_vamb_gen_metadata,
    binning_vamb_prepare_all,
    binning_vamb_postprocess,
    binning_vamb_all