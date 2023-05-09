MULTIBINING_INDEX = {}
for binning_group in SAMPLES_BINNING_GROUP_LIST:
    MULTIBINING_INDEX[binning_group] = {}
    assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, binning_group))
    count = 0
    for assembly_group in assembly_groups:
        count += 1
        MULTIBINING_INDEX[binning_group][assembly_group] = f'''S{count}'''


# reference: https://github.com/RasmussenLab/vamb/blob/master/workflow/vamb.snake.conda.py
rule binning_vamb_combine_scaftigs:
    input:
        lambda wildcards: metapi.get_samples_scaftigs(wildcards, SAMPLES, config["output"]["assembly"])
    output:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_combine_scaftigs/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_combine_scaftigs/{binning_group}.{assembler}.txt")
    params:
        min_contig = config["params"]["binning"]["vamb"]["min_contig"]
    conda:
        config["envs"]["vamb"]
    shell:
        '''
        set +e
        concatenate.py {output} {input} -m {params.min_contig} 2> {log}

        exitcode=$?
        echo "concatenate.py exit code is: $exitcode" >>{log} 2>&1

        gzip -tv {output} 2>&1 | grep -oEi "OK"
        grepcode=$?
        if [ $grepcode -eq 0 ];
        then
            exit 0
        else
            rm -rf {output}
            echo "{output} is not completed." >>{log} 2>&1
            exit $grepcode
        fi
        '''


rule binning_vamb_gen_metadata:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
    output:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz")
    params:
        binning_group = "{binning_group}"
    run:
        import sys
        import gzip

        assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))

        with gzip.open(output[0], 'wt') as oh:
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


localrules:
    binning_vamb_gen_metadata


rule binning_vamb_dict_scaftigs:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
    output:
        os.path.join(
            config["output"]["alignment"],
            "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.dict")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_dict_scaftigs/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_dict_scaftigs/{binning_group}.{assembler}.txt")
    conda:
        config["envs"]["align"]
    shell:
        '''
        samtools dict {input} | \
        cut -f1-3 \
        >{output} \
        2> {log}
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
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_index_scaftigs/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_index_scaftigs/{binning_group}.{assembler}.txt")
    params:
        index_size = config["params"]["binning"]["vamb"]["index_size"]
    threads:
        config["params"]["alignment"]["threads"]
    conda:
        config["envs"]["align"]
    shell:
        '''
        minimap2 -t {threads} \
        -I {params.index_size} \
        -d {output} {input} \
        >{log} 2>&1
        '''


rule binning_vamb_align_scaftigs:
    input:
        reads = os.path.join(SAMPLESDIR, "reads/{sample}/{sample}.json"),
        scaftigs_index = os.path.join(
            config["output"]["alignment"],
            "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.minimap2.mmi"),
        scaftigs_dict = os.path.join(
            config["output"]["alignment"],
            "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.dict")
    output:
        stats = os.path.join(
            config["output"]["alignment"],
            "report/flagstat_minimap2/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.flagstat"),
        bam = temp(os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam")),
        bai = temp(os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam.bai"))
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_align_scaftigs/{binning_group}.{assembler}.{sample}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_align_scaftigs/{binning_group}.{assembler}.{sample}.txt")
    priority:
        28
    threads:
        config["params"]["alignment"]["threads"]
    conda:
        config["envs"]["align"]
    shell:
        '''
        OUTDIR=$(dirname {output.bam})
        OUTPE=$OUTDIR/pe
        OUTSE=$OUTDIR/se

        STATSDIR=$(dirname {output.stats})
        rm -rf $STATSDIR
        mkdir -p $STATSDIR

        rm -rf $OUTDIR
        mkdir -p $OUTDIR

        R1=$(jq -r -M '.PE_FORWARD' {input.reads} | sed 's/^null$//g')
        R2=$(jq -r -M '.PE_REVERSE' {input.reads} | sed 's/^null$//g')
        RS=$(jq -r -M '.SE' {input.reads} | sed 's/^null$//g')

        if [ "$R1" != "" ];
        then
            mkdir -p $OUTPE
            STATSPE=$STATSDIR/minimap2.pe.flagstat

            minimap2 \
            -t {threads} \
            -ax sr \
            {input.scaftigs_index} \
            $R1 $R2 \
            -N 5 2> {log} | \
            tee >(samtools flagstat \
            -@4 - > $STATSPE) | \
            grep -v "^@" | \
            cat {input.scaftigs_dict} - | \
            samtools view -F 3584 -b - | \
            samtools sort \
            -m 3G -@4 \
            -T $OUTPE/temp \
            -O BAM \
            -o $OUTPE/sorted.bam 2>> {log}
        fi

        if [ "$RS" != "" ];
        then
            mkdir -p $OUTSE
            STATSSE=$STATSDIR/minimap2.se.flagstat

            minimap2 \
            -t {threads} \
            -ax sr \
            {input.scaftigs_index} \
            $RS \
            -N 5 2>> {log} | \
            tee >(samtools flagstat \
            -@4 - > $STATSSE) | \
            grep -v "^@" | \
            cat {input.scaftigs_dict} - | \
            samtools view -F 3584 -b - | \
            samtools sort \
            -m 3G -@4 \
            -T $OUTSE/temp \
            -O BAM \
            -o $OUTSE/sorted.bam 2>> {log}
        fi

        if [ -s $OUTPE/sorted.bam ] && [ -s $OUTSE/sorted.bam ];
        then
            samtools merge \
            -l 6 \
            -O BAM -o {output.bam} \
            $OUTPE/sorted.bam \
            $OUTSE/sorted.bam \
            2>>{log}

            rm -rf $OUTPE $OUTSE

            samtools flagstat \
            -@{threads} {output.bam} \
            > {output.stats} \
            2>>{log}

        elif [ -s $OUTPE/sorted.bam ];
        then
            mv $OUTPE/sorted.bam {output.bam}
            mv $STATSPE {output.stats}
            rm -rf $OUTPE

        elif [ -s $OUTSE/sorted.bam ];
        then
            mv $OUTSE/sorted.bam {output.bam}
            mv $STATSSE {output.stats}
            rm -rf $OUTSE
        fi

        samtools index \
        -@{threads} \
        {output.bam} {output.bai} \
        2>> {log}
        '''


rule binning_vamb_align_scaftigs_report:
    input:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat_minimap2/{binning_group}.{{assembler}}/{sample}/{sample}.align2merged_scaftigs.flagstat"),
            zip,
            binning_group=ALIGNMENT_GROUP["binning_group"],
            sample=ALIGNMENT_GROUP["sample_id"])
    output:
        flagstat = os.path.join(config["output"]["alignment"], "report/alignment_flagstat_{assembler}_minimap2.tsv")
    run:
        input_list = [str(i) for i in input]
        output_str = str(output)
        metapi.flagstats_summary(input_list, 2, output=output.flagstat)


rule binning_vamb_coverage:
    input:
        bam = os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam")
    output:
        jgi = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.jgi.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_coverage/{binning_group}.{assembler}.{sample}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_coverage/{binning_group}.{assembler}.{sample}.txt")
    priority:
        30
    conda:
        config["envs"]["metabat2"]
    threads:
        config["params"]["binning"]["threads"]
    shell:
        '''
        JGI={output.jgi}

        jgi_summarize_bam_contig_depths \
        --noIntraDepthVariance \
        --outputDepth ${{JGI%.gz}} \
        {input.bam} \
        >{log} 2>&1

        pigz -f -p {threads} ${{JGI%.gz}}
        '''


rule binning_vamb_gen_abundance_matrix:
    input:
        jgi = lambda wildcards: expand(os.path.join(
            config["output"]["binning"],
            "coverage/{{binning_group}}.{{assembler}}/{sample}.align2merged_scaftigs.jgi.gz"),
            sample=sorted(metapi.get_samples_id_by_binning_group(SAMPLES, wildcards.binning_group)))
    output:
        matrix = os.path.join(
            config["output"]["binning"],
            "matrix/{binning_group}.{assembler}.align2merged_scaftigs.jgi.abundance.matrix.tsv.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb_gen_abundance_matrix/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_gen_abundance_matrix/{binning_group}.{assembler}.txt")
    run:
        metapi.combine_jgi(input.jgi, output.matrix)


rule binning_vamb_prepare_all:
    input:
        expand([
            os.path.join(
                config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
            os.path.join(config["output"]["assembly"],
                "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz"),
            os.path.join(
                config["output"]["binning"],
                "matrix/{binning_group}.{assembler}.align2merged_scaftigs.jgi.abundance.matrix.tsv.gz"),
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
            "matrix/{binning_group}.{assembler}.align2merged_scaftigs.jgi.abundance.matrix.tsv.gz")
    output:
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_vamb/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb/{binning_group}.{assembler}.log")
    params:
        outdir = os.path.join(config["output"]["binning"], "mags_vamb/{binning_group}.{assembler}"),
        min_contig = config["params"]["binning"]["vamb"]["min_contig"],
        min_fasta = config["params"]["binning"]["vamb"]["min_fasta"],
        cuda = "--cuda" if config["params"]["binning"]["vamb"]["cuda"] else "",
        cuda_module = config["params"]["binning"]["vamb"]["cuda_module"],
        use_cuda_module = int(config["params"]["binning"]["vamb"]["use_cuda_module"]),
        allow_small_scaftigs = 1 if config["params"]["binning"]["vamb"]["allow_small_scaftigs"] else 0,
        external_params = config["params"]["binning"]["vamb"]["external_params"]
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

        MATRIX={input.matrix}
        pigz -dkf $MATRIX

        vamb \
        {params.cuda} \
        -p {threads} \
        --outdir {params.outdir} \
        --fasta {input.scaftigs} \
        --jgi ${{MATRIX%.gz}} \
        -o C \
        -m {params.min_contig} \
        --minfasta {params.min_fasta} \
        {params.external_params} \
        >> {log} 2>&1

        rm -rf ${{MATRIX%.gz}}


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
        metadata = os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz"),
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}/binning_done")
    output:
        metadata = os.path.join(
            config["output"]["binning"],
            "mags_vamb/{binning_group}.{assembler}/bins_{assembly_group}/cluster.metadata.tsv.gz"),
        binning_done = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/vamb/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_postprocess/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_vamb_postprocess/{binning_group}.{assembly_group}.{assembler}.txt")
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
                shell(f'''pigz -f {fna}''')
                bin_index += 1
                # bin_id = os.path.basename(fna).split(".")[0]
                # bin_id = os.path.basename(fna).split(".")[0].split("C")[-1]
                fna_dist = os.path.join(outdir, f'''{params.binning_group}.{params.assembly_group}.{params.assembler}.vamb.bin.{bin_index}.fa.gz''')
                metadata.append((os.path.abspath(fna) + ".gz", os.path.abspath(fna_dist)))
                shell(f'''zcat {fna}.gz | seqkit replace -p "^S\d+C" | pigz -cf > {fna_dist}''')

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
                    "mags_vamb/{binning_group}.{assembler}/bins_{assembly_group}/cluster.metadata.tsv.gz"),
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
    binning_vamb_prepare_all,
    binning_vamb_all
