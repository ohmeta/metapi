rule filter_rename_prediction:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        cds = os.path.join(config["results"]["cobinning"]["cds"],
                           "{sample}/{sample}.{assembler}.cds.fa.gz"),
        gff = os.path.join(config["results"]["cobinning"]["cds"],
                           "{sample}/{sample}.{assembler}.cds.gff.gz")
    params:
        cds = os.path.join(config["results"]["cobinning"]["cds"],
                           "{sample}/{sample}.{assembler}.cds.fa"),
        gff = os.path.join(config["results"]["cobinning"]["cds"],
                           "{sample}/{sample}.{assembler}.cds.gff"),
        length = config["params"]["cobinning"]["scaftigs_length"],
        id = "{sample}"
    threads:
        config["params"]["cobinning"]["threads"]
    shell:
        '''
        seqkit seq --min-len {params.length} --threads {threads} --quiet {input} | 
        seqkit replace --pattern "^(.+)" --replacement '{params.id}_$1' |
        prodigal -d {params.cds} -o {params.gff} -f gff -p meta -q
        pigz {params.cds}
        pigz {params.gff}
        '''


rule vsearch_clust_cds:
    input:
        cds = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.fa.gz")
    output:
        uc = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.uc.gz")
    params:
        uc = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.uc"),
        identity = config["params"]["cobinning"]["vsearch"]["identity"]
    log:
        os.path.join(config["logs"]["cobinning"], "{sample}.{assembler}.vsearch_clust.cds.log")
    threads:
        config["params"]["cobinning"]["threads"]
    shell:
        '''
        pigz -d -c {input.cds} |
        vsearch --threads {threads} --cluster_fast - --id {params.identity} --uc {params.uc} 2> {log}
        pigz {params.uc}
        '''


rule choose_cds_marker:
    input:
        cds = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.fa.gz"),
        uc = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.uc.gz")
    output:
        cds = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.marker.fa.gz")
    log:
        os.path.join(config["logs"]["cobinning"], "{sample}.{assembler}.cds.marker.log")
    run:
        from Bio import SeqIO
        import gzip

        logh = open(log, 'w')

        clust = set()
        hit = set()
        with gzip.open(input.uc, 'r') as uch:
            for line in uch:
                item = line.strip().split('\t')
                if 'H' in item[0]:
                    hit.add(item[-1])
                    hit.add(item[-2])
                elif 'C' in item[0]:
                    clust.add(item[-2])
        # understand it
        marker_cds = clust - hit

        with gzip.open(input.cds, 'rt') as ih, gzip.open(output.cds, 'wb') as oh:
            records = [r for r in SeqIO.parse(ih, "fasta") if r.id in marker_cds]
            cds_count = SeqIO.write(records, oh, "fasta")
            logh.write("Saved %i marker genes to %s\n" % (cds_count, output.cds))
        logh.close()


rule index_marker_cds:
    input:
        os.path.join(config["results"]["cobinning"]["cds"],
                     "{sample}/{sample}.{assembler}.cds.marker.fa.gz")
    output:
        expand("{prefix}.{suffix}",
               prefix=os.path.join(config["results"]["cobinning"]["cds"],
                                   "{{sample}}/{{sample}}.{{assembler}}.cds.marker.fa.gz"),
               suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    log:
        os.path.join(config["logs"]["cobinning"], "{sample}.{assembler}.cds.marker.index.log")
    params:
        prefix = os.path.join(config["results"]["cobinning"]["cds"],
                              "{sample}/{sample}.{assembler}.cds.marker.fa.gz")
    shell:
        '''
        bowtie2-build {input} {params.prefix} >{log} 2>&1
        '''


rule alignment_to_marker_cds:
    input:
        reads = clean_reads,
        db = expand("{prefix}.{suffix}",
                    prefix=os.path.join(config["results"]["cobinning"]["cds"],
                                        "{{sample_}}/{{sample_}}.{{assembler}}.cds.marker.fa.gz"),
                    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    output:
        os.path.join(config["results"]["cobinning"]["depth"],
                     "{sample}/{sample}.{sample_}.{assembler}.metabat2.depth.txt")
    params:
        index = os.path.join(config["results"]["cobinning"]["cds"], "{sample_}/{sample_}.{assembler}.cds.marker.fa.gz"),
        bam = os.path.join(config["results"]["cobinning"]["depth"],
                           "{sample}/{sample}.{sample_}.{assembler}.sorted.bam")
    threads:
        config["params"]["cobinning"]["threads"]
    shell:
         '''
         bowtie2 --threads {threads} -x {params.index} -1 {input.reads[0]} -2 {input.reads[1]} |
         samtools sort -@{threads} -O BAM -o {params.bam} -
         jgi_summarize_bam_contig_depths --outputDepth {output} {params.bam}
         rm -rf {params.bam}
         '''

