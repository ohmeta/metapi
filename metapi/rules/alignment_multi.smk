rule alignment_scaftigs_combined:
    input:
        scaftigs = lambda wildcards: metapi.get_samples_scaftigs(wildcards, SAMPLES, config["output"]["assembly"])
    output:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/alignment_scaftigs_combined/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/alignment_scaftigs_combined/{binning_group}.{assembler}.txt")
    params:
        min_contig = config["params"]["binning"]["min_contig_len_bp"],
        separator = config["params"]["binning"]["separator"]
    threads:
        config["params"]["alignment"]["threads"]
    conda:
        config["envs"]["report"]
    shell:
        '''
        rm -rf {output.scaftigs}

        SCAFTIGSGZ={output.scaftigs}
        SCAFTIGS=${{SCAFTIGSGZ%.gz}}

        for FASTA in {input.scaftigs}
        do
            ASMINDEX=$(basename $FASTA | sed 's#.scaftigs.fa.gz##g')   

            bioawk \
            -v scaftigsid=$ASMINDEX \
            -v separator={params.separator} \
            -v mincontig={params.min_contig} \
            -c fastx '{{if(length($seq) >= mincontig){{print ">" scaftigsid separator $name;print $seq}}}}' $FASTA \
            >> $SCAFTIGS \
            2>> {log}
        done
        
        pigz -p {threads} $SCAFTIGS 
        '''
    # too slow
    #run:
    #    import os
    #    import sys
    #    import gzip
    #    from Bio import SeqIO

    #    scaftigs_sorted = sorted(input.scaftigs) 

    #    with gzip.open(output.scaftigs, "wt") as oh:
    #        for scaftigs in scaftigs_sorted:
    #            sample_name = os.path.basename(scaftigs).replace(".scaftigs.fa.gz", "")
    #            with gzip.open(scaftigs, "rt") as ih:
    #                for rc in SeqIO.parse(ih, "fasta"):
    #                    if len(rc) >= params.min_contig:
    #                        seq_id = rc.id
    #                        rc.id = f'{sample_name}{params.separator}{seq_id}'
    #                        SeqIO.write(rc, oh, "fasta")
   

rule alignment_scaftigs_combined_dict:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz")
    output:
        scaftigs_dict = os.path.join(
            config["output"]["alignment"],
            "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.dict"),
        scaftigs_headers = os.path.join(
            config["output"]["alignment"],
            "index_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.headers.txt"),
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/alignment_scaftigs_combined_dict/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/alignment_scaftigs_combined_dict/{binning_group}.{assembler}.txt")
    conda:
        config["envs"]["align"]
    shell:
        '''
        samtools dict {input} | \
        cut -f1-3 \
        > {output.scaftigs_dict} \
        2> {log}

        tail -n +2 {output.scaftigs_dict} | \
        cut -f2,3 \
        > {output.scaftigs_headers} \
        2>> {log}
        '''


rule alignment_scaftigs_combined_index:
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
            "logs/alignment_scaftigs_combined_index/{binning_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/alignment_scaftigs_combined_index/{binning_group}.{assembler}.txt")
    params:
        index_size = config["params"]["binning"]["combined_index_size"]
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


rule alignment_scaftigs_reads_multi:
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
            "report/flagstat_minimap2/{binning_group}.{assembler}/{sample}.align2merged_scaftigs.flagstat"),
        bam = os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam"),
        bai = os.path.join(
            config["output"]["alignment"],
            "bam_merged/{binning_group}.{assembler}/{sample}/{sample}.align2merged_scaftigs.sorted.bam.bai")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/alignment_scaftigs_reads_multi/{binning_group}.{assembler}.{sample}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/alignment_scaftigs_reads_multi/{binning_group}.{assembler}.{sample}.txt")
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
            -o $OUTPE/sorted.bam
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
            -o $OUTSE/sorted.bam
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


rule alignment_scaftigs_reads_multi_report:
    input:
        expand(
            os.path.join(
                config["output"]["alignment"],
                "report/flagstat_minimap2/{binning_group}.{{assembler}}/{sample}.align2merged_scaftigs.flagstat"),
            zip,
            binning_group=ALIGNMENT_GROUP["binning_group"],
            sample=ALIGNMENT_GROUP["sample_id"])
    output:
        flagstat = os.path.join(config["output"]["alignment"], "report/alignment_flagstat_{assembler}_minimap2.tsv")
    run:
        input_list = [str(i) for i in input]
        output_str = str(output)
        metapi.flagstats_summary(input_list, 2, output=output.flagstat)


rule alignment_scaftigs_reads_multi_report_all:
    input:
        expand(os.path.join(config["output"]["alignment"], "report/alignment_flagstat_{assembler}_minimap2.tsv"),
        assembler=ASSEMBLERS)


rule alignment_all:
    input:
        rules.alignment_base_depth_all.input,
        rules.alignment_scaftigs_reads_report_all.input,
        rules.alignment_scaftigs_reads_multi_report_all.input


localrules:
    alignment_scaftigs_reads_all,
    alignment_base_depth_all,
    alignment_scaftigs_reads_report_all,
    alignment_scaftigs_reads_multi_report_all,
    alignment_all
