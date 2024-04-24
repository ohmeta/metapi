MULTIBINING_INDEX = {}
for binning_group in SAMPLES_BINNING_GROUP_LIST:
    MULTIBINING_INDEX[binning_group] = {}
    assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, binning_group))
    count = 0
    for assembly_group in assembly_groups:
        count += 1
        MULTIBINING_INDEX[binning_group][assembly_group] = f'''S{count}'''


rule binning_vamb_gen_abundance_mask:
    input:
        headers = os.path.join(
            config["output"]["alignment"],
            "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.headers.txt")
    output:
        mask_refhash = os.path.join(
            config["output"]["binning"],
            "matrix/{binning_group}.{assembler}.mask_refhash.npz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_gen_abundance_mask/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_gen_abundance_mask/{binning_group}.{assembler}.txt")
    params:
        script = os.path.join(WRAPPER_DIR, "vamb", "abundances_mask.py"),
        min_contig = config["params"]["binning"]["min_contig_len_bp"]
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["vamb"]
    shell:
        '''
        python {params.script} \
        --h {input.headers} \
        --msk {output.mask_refhash} \
        --minsize {params.min_contig} \
        2> {log}
        '''


rule binning_vamb_gen_abundance_samples:
    input:
        bam = os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam"),
        bai = os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam.bai"),
        mask_refhash = os.path.join(
            config["output"]["binning"],
            "matrix/{binning_group}.{assembler}.mask_refhash.npz")
    output:
        abundance = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.npz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_gen_abundance_samples/{binning_group}.{assembler}.{sample}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_gen_abundance_samples/{binning_group}.{assembler}.{sample}.txt")
    params:
        script = os.path.join(WRAPPER_DIR, "vamb", "write_abundances.py"),
        min_identity = config["params"]["binning"]["vamb"]["min_identity"]
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["vamb"]
    shell:
        '''
        python {params.script} \
        --msk {input.mask_refhash} \
        --b {input.bam} \
        --min_id {params.min_identity} \
        --out {output.abundance} \
        2> {log}
        '''


rule binning_vamb_gen_abundance_matrix:
    input:
        abundances = lambda wildcards: expand(os.path.join(
            config["output"]["binning"],
            "coverage/{{binning_group}}.{{assembler}}/{sample}.align2merged_scaftigs.npz"),
            sample=sorted(metapi.get_samples_id_by_binning_group(SAMPLES, wildcards.binning_group))),
        mask_refhash = os.path.join(
            config["output"]["binning"],
            "matrix/{binning_group}.{assembler}.mask_refhash.npz")
    output:
        matrix = os.path.join(
            config["output"]["binning"],
            "matrix/{binning_group}.{assembler}.abundance.matrix.npz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_gen_abundance_matrix/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_gen_abundance_matrix/{binning_group}.{assembler}.txt")
    params:
        script = os.path.join(WRAPPER_DIR, "vamb", "create_abundances.py"),
        min_identity = config["params"]["binning"]["vamb"]["min_identity"]
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["vamb"]
    shell:
        '''
        python {params.script} \
        --msk {input.mask_refhash} \
        --ab {input.abundances} \
        --min_id {params.min_identity} \
        --out {output} \
        2> {log}
        '''


rule binning_vamb_prepare_all:
    input:
        expand([
            os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
            os.path.join(
                config["output"]["binning"],
                "matrix/{binning_group}.{assembler}.abundance.matrix.npz"),
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
            "matrix/{binning_group}.{assembler}.abundance.matrix.npz")
    output:
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}.{vamber}/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_run_{vamber}/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_run_{vamber}/{binning_group}.{assembler}.log")
    wildcard_constraints:
        vamber="[a]?vamb"
    params:
        outdir = os.path.join(config["output"]["binning"], "mags_vamb/{binning_group}.{assembler}.{vamber}"),
        min_contig = config["params"]["binning"]["min_contig_len_bp"],
        min_fasta = config["params"]["binning"]["min_bin_len_kbp"] * 1000,
        cuda = "--cuda" if config["params"]["binning"]["vamb"]["cuda"] else "",
        cuda_module = config["params"]["binning"]["vamb"]["cuda_module"],
        use_cuda_module = int(config["params"]["binning"]["vamb"]["use_cuda_module"]),
        allow_small_scaftigs = 1 if config["params"]["binning"]["vamb"]["allow_small_scaftigs"] else 0,
        binner = lambda wildcards: "default" if wildcards.vamber == "vamb" else "avamb",
        seed = config["params"]["seed"],
        external_params = config["params"]["binning"]["vamb"]["external_params"],
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["vamb"]
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
            lspci | grep -oEi nvidia >> {log} 2>&1
            grepcode=$?
            if [ $grepcode -ne 0 ];
            then
                echo "No NVIDIA GPU detected, please change vamb::use_cuda to false and rerun the pipeline."
                exit 0
            else
                echo "NVIDIA GPU detected, you specific vamb::use_cuda to true, great!"
                which python >> {log} 2>&1
                which vamb >> {log} 2>&1

                python -c 'import torch;print(torch.__file__)' >> {log} 2>&1
                python -c 'import torch;print(f"Torch CUDA: {{torch.cuda.is_available()}}")' >> {log} 2>&1
                python -c 'from torch.utils.cpp_extension import CUDA_HOME;print(CUDA_HOME)' >> {log} 2>&1
                python -c 'import os; print(os.environ.get("CUDA_PATH"))' >> {log} 2>&1
            fi
        fi

        vamb bin {params.binner} \
        {params.cuda} \
        -p {threads} \
        --seed {params.seed} \
        --outdir {params.outdir} \
        --fasta {input.scaftigs} \
        --rpkm {input.matrix} \
        -o C \
        -m {params.min_contig} \
        --minfasta {params.min_fasta} \
        {params.external_params} \
        >> {log} 2>&1

        if [ -f {params.outdir}/clusters.tsv ];
        then
            echo "Running vamb completed" >> {log} 2>&1
            echo "Touch binning_done" >> {log} 2>&1
            touch {output.binning_done}
            exit 0
        else
            echo "No bins generated, please check {log}"
            echo "No bins generated, please check {log}" >>{log}
            exit 1
        fi
        '''


rule binning_vamb_postprocess:
    input:
        #metadata = os.path.join(
        #    config["output"]["assembly"],
        #    "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz"),
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}.{vamber}/binning_done")
    output:
        metadata = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}.{vamber}/bins_{assembly_group}/cluster.metadata.tsv.gz"),
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/{vamber}/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_postprocess_{vamber}/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_postprocess_{vamber}/{binning_group}.{assembly_group}.{assembler}.txt")
    wildcard_constraints:
        vamber="[a]?vamb"
    params:
        binning_group = "{binning_group}",
        assembly_group = "{assembly_group}",
        assembler = "{assembler}",
        vamber = "{vamber}",
        separator = config["params"]["binning"]["separator"]
    run:
        from glob import glob
        import os
        import sys
        import pandas as pd

        # binning_assembly_metadata = pd.read_csv(input.metadata, sep="\t").set_index("binning_assembly_group")
        # assembly_index = binning_assembly_metadata.loc[f'''{params.binning_group}.{params.assembly_group}''', "vamb_id"]

        #assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))
        #assembly_index = int(assembly_groups.index(params.assembly_group)) + 1
        #assembly_index = f'''S{assembly_index}'''
        ## Double check
        #if assembly_index != MULTIBINING_INDEX[params.binning_group][params.assembly_group]:
        #    sys.exit("assembly_group index error")

        assembly_index = f'''{params.binning_group}.{params.assembly_group}.{params.assembler}'''
        metadata = []

        outdir = os.path.dirname(output.binning_done)
        mags_dir = os.path.dirname(input.binning_done)
        os.makedirs(outdir, exist_ok=True)
        bin_index = 0

        if os.path.exists(f'{mags_dir}/bins'):
            fna_list = sorted(glob(f'{mags_dir}/bins/{assembly_index}{params.separator}*.fna'))

            for fna in fna_list:
                shell(f'''pigz -f {fna}''')
                bin_index += 1
                # bin_id = os.path.basename(fna).split(".")[0]
                # bin_id = os.path.basename(fna).split(".")[0].split("C")[-1]
                fna_dist = os.path.join(outdir, f'''{assembly_index}.{params.vamber}.bin.{bin_index}.fa.gz''')
                metadata.append((os.path.abspath(fna) + ".gz", os.path.abspath(fna_dist)))
                shell(f'''zcat {fna}.gz | seqkit replace -p "^S\\w+{params.separator}" | pigz -cf > {fna_dist}''')

        shell(f'''touch {output.binning_done}''')

        pd.DataFrame(metadata, columns=["vamb_bin", "vamb_postprocess_bin"])\
            .to_csv(output.metadata, sep='\t', index=False)


localrules:
    binning_vamb_postprocess


if config["params"]["binning"]["vamb"]["do"]:
    rule binning_vamb_all:
        input:
            rules.binning_vamb_prepare_all.input,
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "mags_vamb/{binning_group}.{assembler}.{vamber}/{results}"),
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS,
                vamber=BINNERS_VAMB,
                results=[
                    #"clusters.tsv",
                    #"latent.npz",
                    #"lengths.npz",
                    #"log.txt",
                    #"model.pt",
                    #"mask.npz",
                    #"tnf.npz",
                    "binning_done"]),
            expand(expand([
                os.path.join(
                    config["output"]["binning"],
                    "mags_vamb/{binning_group}.{assembler}.{{vamber}}/bins_{assembly_group}/cluster.metadata.tsv.gz"),
                os.path.join(
                    config["output"]["binning"],
                    "mags/{binning_group}.{assembly_group}.{assembler}/{{vamber}}/binning_done")],
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
                vamber=BINNERS_VAMB)

else:
    rule binning_vamb_all:
        input:


rule binning_semibin_multi_prepare_bam:
    input:
        bam = lambda wildcards: metapi.get_samples_bax_multi(wildcards, SAMPLES, config["output"]["alignment"], "bam"),
        bai = lambda wildcards: metapi.get_samples_bax_multi(wildcards, SAMPLES, config["output"]["alignment"], "bam.bai")
    output:
        bam = os.path.join(
            config["output"]["alignment"],
            "bam_merged_by_assembly_group/{binning_group}.{assembly_group}.{assembler}/{assembly_group}.sorted.bam"),
        bai = os.path.join(
            config["output"]["alignment"],
            "bam_merged_by_assembly_group/{binning_group}.{assembly_group}.{assembler}/{assembly_group}.sorted.bam.bai")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_semibin_multi_prepare_bam/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_semibin_multi_prepare_bam/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        bam_dir = os.path.join(
            config["output"]["alignment"],
            "bam_merged_by_assembly_group/{binning_group}.{assembly_group}.{assembler}")
    conda:
        config["envs"]["align"]
    threads:
        config["params"]["alignment"]["threads"]
    shell:
        '''
        rm -rf {params.bam_dir}
        mkdir -p {params.bam_dir}

        inputarray=({input.bam})
        inputlen=${{#inputarray[@]}}

        if [ $inputlen -eq 1 ];
        then 
            BAM=$(realpath {input.bam[0]})
            BAI=$(realpath {input.bai[0]})
            ln -s $BAM {output.bam}
            ln -s $BAI {output.bai}
        else
            samtools merge \
            -l 6 \
            -O BAM -o {output.bam}.merged {input.bam} \
            2> {log}

            samtools sort \
            -m 3G \
            -@4 \
            -T {output.bam}.temp \
            -O BAM -o {output.bam} {output.bam}.merged \
            2>> {log}

            rm -rf {output.bam}.merged
            rm -rf {output.bam}.temp*

            samtools index \
            -@{threads} \
            {output.bam} {output.bai} \
            2>> {log}
        fi
        '''
 

rule binning_semibin_multi_generate_sequence_features:
    input:
        scaftigs = os.path.join(config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
        bam = lambda wildcards: expand(os.path.join(
            config["output"]["alignment"],
            "bam_merged_by_assembly_group/{{binning_group}}.{assembly_group}.{{assembler}}/{assembly_group}.sorted.bam"),
            assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group))),
        bai = lambda wildcards: expand(os.path.join(
            config["output"]["alignment"],
            "bam_merged_by_assembly_group/{{binning_group}}.{assembly_group}.{{assembler}}/{assembly_group}.sorted.bam.bai"),
            assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
    output:
        #data = expand(os.path.join(
        #    config["output"]["binning"],
        #    "mags_semibin_multi/{{binning_group}}.{{assembler}}.semibin_multi/samples/{assembly_group}/data.csv"),
        #    assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, "{binning_group}"))),
        #data_split = expand(os.path.join(
        #    config["output"]["binning"],
        #    "mags_semibin_multi/{{binning_group}}.{{assembler}}.semibin_multi/samples/{assembly_group}/data_split.csv"),
        #    assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, "{binning_group}")))
        done = os.path.join(
            config["output"]["binning"],
            "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/binning_generate_sequence_features_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_semibin_multi_generate_sequence_features/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_semibin_multi_generate_sequence_features/{binning_group}.{assembler}.txt")
    params:
        separator = config["params"]["binning"]["separator"],
        min_len = config["params"]["binning"]["min_contig_len_bp"],
        out_dir = os.path.join(
            config["output"]["binning"],
            "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi")
    conda:
        config["envs"]["semibin"]
    threads:
        config["params"]["binning"]["threads"]
    shell:
        '''
        SemiBin2 generate_sequence_features_multi \
        --input-fasta {input.scaftigs} \
        --input-bam {input.bam} \
        --compression gz \
        --separator {params.separator} \
        --min-len {params.min_len} \
        --threads {threads} \
        --output {params.out_dir} \
        >{log} 2>&1
        '''


rule binning_semibin_multi_generate_sequence_features_all:
    input:
        expand(os.path.join(
            config["output"]["binning"],
            "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/binning_generate_sequence_features_done"),
            binning_group=SAMPLES_BINNING_GROUP_LIST,
            assembler=ASSEMBLERS)


# it seems that it is same with binning_semibin_single_generate_cannot_links
rule binning_semibin_multi_generate_cannot_links:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/cannot/cannot.txt")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_semibin_multi_generate_cannot_links/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_semibin_multi_generate_cannot_links/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        min_len = config["params"]["binning"]["min_contig_len_bp"],
        reference_db = config["params"]["binning"]["semibin"]["reference_db"],
        out_dir = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi")
    conda:
        config["envs"]["semibin"]
    threads:
        config["params"]["binning"]["threads"]
    shell:
        '''
        SemiBin2 generate_cannot_links \
        --input-fasta {input.scaftigs} \
        --compression gz \
        --reference_db {params.reference_db} \
        --min-len {params.min_len} \
        --threads {threads} \
        --output {params.out_dir} \
        >{log} 2>&1
        '''


if  config["params"]["binning"]["semibin"]["train_mode"] == "self":
    rule binning_semibin_multi_train:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
            #data = os.path.join(
            #    config["output"]["binning"],
            #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
            #data_split = os.path.join(
            #    config["output"]["binning"],
            #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv")
            done = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/binning_generate_sequence_features_done")
        output:
            #model = os.path.join(
            #    config["output"]["binning"],
            #    "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/model.h5")
            done = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/train_done")
        log:
            os.path.join(
                config["output"]["binning"],
                "logs/binning_semibin_multi_train/{binning_group}.{assembly_group}.{assembler}.log")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/binning_semibin_multi_train/{binning_group}.{assembly_group}.{assembler}.txt")
        params:
            min_len = config["params"]["binning"]["min_contig_len_bp"],
            train_mode = config["params"]["binning"]["semibin"]["train_mode"],
            train_batch_size = config["params"]["binning"]["semibin"]["train_batch_size"],
            train_epochs = config["params"]["binning"]["semibin"]["train_epochs"],
            random_seed = config["params"]["seed"],
            engine = config["params"]["binning"]["semibin"]["engine"],
            orf_finder = config["params"]["binning"]["semibin"]["orf_finder"],
            out_dir = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi"),
            data = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
            data_split = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv")
        conda:
            config["envs"]["semibin"]
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            SemiBin2 train_self \
            --data {params.data} \
            --data-split {params.data_split} \
            --batch-size {params.train_batch_size} \
            --epochs {params.train_epochs} \
            --random-seed {params.random_seed} \
            --engien {params.engine} \
            --output {params.out_dir} \
            --threads {threads} \
            >{log} 2>&1

            touch {output.done}
            '''


elif  config["params"]["binning"]["semibin"]["train_mode"] == "semi":
    rule binning_semibin_multi_train:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
            #data = os.path.join(
            #    config["output"]["binning"],
            #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
            #data_split = os.path.join(
            #    config["output"]["binning"],
            #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv"),
            cannot = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/cannot/cannot.txt"),
            done = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/binning_generate_sequence_features_done")
        output:
            #model = os.path.join(
            #    config["output"]["binning"],
            #    "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/model.h5")
            done = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/train_done")
        log:
            os.path.join(
                config["output"]["binning"],
                "logs/binning_semibin_multi_train/{binning_group}.{assembly_group}.{assembler}.log")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/binning_semibin_multi_train/{binning_group}.{assembly_group}.{assembler}.txt")
        params:
            min_len = config["params"]["binning"]["min_contig_len_bp"],
            train_mode = config["params"]["binning"]["semibin"]["train_mode"],
            train_batch_size = config["params"]["binning"]["semibin"]["train_batch_size"],
            train_epochs = config["params"]["binning"]["semibin"]["train_epochs"],
            random_seed = config["params"]["seed"],
            engine = config["params"]["binning"]["semibin"]["engine"],
            orf_finder = config["params"]["binning"]["semibin"]["orf_finder"],
            out_dir = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi"),
            data = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
            data_split = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv")
        conda:
            config["envs"]["semibin"]
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            SemiBin2 train_semi \
            --input-fasta {input.scaftigs} \
            --cannot-link {input.cannot} \
            --data {params.data} \
            --data-split {params.data_split} \
            --batch-size {params.train_batch_size} \
            --epochs {params.train_epochs} \
            --random-seed {params.random_seed} \
            --engien {params.engine} \
            --orf-finder {params.orf_finder} \
            --output {params.out_dir} \
            --threads {threads} \
            >{log} 2>&1

            touch {output.done}
            '''

else:
    rule binning_semibin_multi_train:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
            #data = os.path.join(
            #    config["output"]["binning"],
            #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
            #data_split = os.path.join(
            #    config["output"]["binning"],
            #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv")
            done = os.path.join(
                config["output"]["binning"],
                "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/binning_generate_sequence_features_done")
        output:
            done = os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/train_done")
        log:
            os.path.join(
                config["output"]["binning"],
                "logs/binning_semibin_multi_train/{binning_group}.{assembly_group}.{assembler}.log")
        benchmark:
            os.path.join(
                config["output"]["binning"],
                "benchmark/binning_semibin_multi_train/{binning_group}.{assembly_group}.{assembler}.txt")
        conda:
            config["envs"]["semibin"]
        threads:
            1
        shell:
            '''
            touch {output.done}
            '''

 
rule binning_semibin_multi_bin:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
        #data = os.path.join(
        #    config["output"]["binning"],
        #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
        #data_split = os.path.join(
        #    config["output"]["binning"],
        #    "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv"),
        done_data = os.path.join(
            config["output"]["binning"],
            "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/binning_generate_sequence_features_done"),
        done_train = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/train_done")
    output:
        os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_semibin_multi_bin/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_semibin_multi_bin/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        min_len = config["params"]["binning"]["min_contig_len_bp"],
        min_fasta = config["params"]["binning"]["min_bin_len_kbp"],
        out_dir = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi"),
        train_mode = config["params"]["binning"]["semibin"]["train_mode"],
        environment = config["params"]["binning"]["semibin"]["environment"],
        seed = config["params"]["seed"],
        engine = config["params"]["binning"]["semibin"]["engine"],
        data = os.path.join(
            config["output"]["binning"],
            "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data.csv"),
        data_split = os.path.join(
            config["output"]["binning"],
            "mags_semibin_multi/{binning_group}.{assembler}.semibin_multi/samples/{assembly_group}/data_split.csv")
    conda:
        config["envs"]["semibin"]
    threads:
        config["params"]["binning"]["threads"]
    shell:
        '''
        if [ "{params.train_mode}" == "self" ] || [ "{params.train_mode}" == "semi" ];
        then
            SemiBin2 bin_short \
            --threads {threads} \
            --input-fasta {input.scaftigs} \
            --data {params.data} \
            --compression gz \
            --model {params.out_dir}/model.h5 \
            --engine {params.engine} \
            --random-seed {params.seed} \
            --min-len {params.min_len} \
            --minfasta-kbs {params.min_fasta} \
            --output {params.out_dir} \
            >{log} 2>&1

            touch {output}

        else
            SemiBin2 bin_short \
            --threads {threads} \
            --input-fasta {input.scaftigs} \
            --data {params.data} \
            --compression gz \
            --environment {params.environment} \
            --engine {params.engine} \
            --random-seed {params.seed} \
            --min-len {params.min_len} \
            --minfasta-kbs {params.min_fasta} \
            --output {params.out_dir} \
            >{log} 2>&1

            touch {output}
        fi
        '''


if config["params"]["binning"]["semibin"]["do"]:
    if "multi" in config["params"]["binning"]["semibin"]["mode"]:
        rule binning_semibin_multi_bin_all:
            input:
                expand(os.path.join(
                    config["output"]["binning"],
                    "mags/{binning_group}.{assembly_group}.{assembler}/semibin_multi/binning_done"),
                    zip,
                    binning_group=ASSEMBLY_GROUPS["binning_group"],
                    assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                    assembler=ASSEMBLY_GROUPS["assembler"])

    else:
        rule binning_semibin_multi_bin_all:
            input:
else:
    rule binning_semibin_multi_bin_all:
        input:


rule binning_semibin_all:
    input:
        rules.binning_semibin_single_bin_all.input,
        rules.binning_semibin_multi_bin_all.input


localrules:
    binning_vamb_prepare_all,
    binning_vamb_all,
    binning_semibin_multi_generate_sequence_features_all,
    binning_semibin_multi_bin_all,
    binning_semibin_all

