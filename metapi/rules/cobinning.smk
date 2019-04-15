

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
        id = lambda wildcards: renamed_id(_samples, wildcards) \
                               if config["params"]["cobinning"]["rename"] else "{sample}",
        assembler = "{assembler}"
    threads:
        config["params"]["cobinning"]["threads"]
    shell:
        '''
        if [[ {params.assembler} == "metaspades" ]]; then
            seqkit seq --min-len {params.length} --threads {threads} --quiet {input} | 
            seqkit replace --threads {threads} --pattern "^\w+?_(\d+)?_.*" --replacement '{params.id}_$1' |
            prodigal -d {params.cds} -o {params.gff} -f gff -p meta -q
            pigz {params.cds}
            pigz {params.gff}
        fi
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

        logh = open(str(log), 'w')

        clust = set()
        hit = set()
        with gzip.open(input.uc, 'rt') as uch:
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


def clean_reads_(wildcards):
    if config["params"]["begin"] == "assembly":
        if config["params"]["type"] == "fastq":
            r1 = get_sample_id_(_samples, wildcards, "fq1")
            r2 = get_sample_id_(_samples, wildcards, "fq2")
            return [r1, r2]
        elif config["params"]["type"] == "sra":
            return expand("{sra2fq}/{sample_}.{read}.fq.gz",
                          sra2fq=config["results"]["sra2fq"],
                          sample=wildcards.sample_,
                          read=[1, 2])
    elif config["params"]["rmhost"]["do"]:
        return expand("{rmhost}/{sample_}.rmhost.{read}.fq.gz",
                      rmhost=config["results"]["rmhost"],
                      sample_=wildcards.sample_,
                      read=["1", "2"])
    else:
        return expand("{trimming}/{sample_}.trimmed.{read}.fq.gz",
                      trimming=config["results"]["trimming"],
                      sample_=wildcards.sample_,
                      read=["1", "2"])


rule get_marker_contigs_depth:
    input:
        reads = clean_reads_,
        db = expand("{prefix}.{suffix}",
                    prefix=os.path.join(config["results"]["cobinning"]["cds"],
                                        "{{sample}}/{{sample}}.{{assembler}}.cds.marker.fa.gz"),
                    suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    output:
        os.path.join(config["results"]["cobinning"]["depth"],
                     "{sample_}/{sample_}.{sample}.{assembler}.metabat2.depth.txt.gz")
    params:
        index = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.marker.fa.gz"),
        depth = os.path.join(config["results"]["cobinning"]["depth"],
                             "{sample_}/{sample_}.{sample}.{assembler}.metabat2.depth.txt")
    threads:
        config["params"]["cobinning"]["threads"]
    shell:
         '''
         bowtie2 --threads {threads} -x {params.index} -1 {input.reads[0]} -2 {input.reads[1]} |
         samtools sort -@{threads} -O BAM - |
         jgi_summarize_bam_contig_depths --outputDepth {params.depth} - 
         pigz {params.depth}
         '''
