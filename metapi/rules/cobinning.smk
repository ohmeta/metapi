def cobinning_input(wildcards):
    if RMHOST_DO:
        return get_reads_(wildcards, "rmhost", False)
    elif TRIMMING_DO:
        return get_reads_(wildcards, "trimming", False)
    else:
        return get_reads_(wildcards, "raw", False)


if config["params"]["cobinning"]["do"]:
    rule cobinning_vsearch_clust_cds:
        input:
            cds = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.ffn")
        output:
            uc = os.path.join(config["output"]["cobinning"],
                              "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.uc.gz")
        log:
            os.path.join(config["output"]["cobinning"],
                         "logs/clust/{sample}.{assembler}.vsearch_clust.cds.log")
        params:
            uc = os.path.join(config["output"]["cobinning"],
                              "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.uc"),
            identity = config["params"]["cobinning"]["vsearch"]["identity"]
        threads:
            config["params"]["cobinning"]["threads"]
        shell:
            '''
            vsearch \
            --threads {threads} \
            --cluster_fast {input.cds} \
            --id {params.identity} \
            --uc {params.uc} \
            2> {log}

            pigz {params.uc}
            '''


    rule cobinning_choose_cds_marker:
        input:
            cds = os.path.join(
                config["output"]["predict"],
                "scaftigs_gene/{sample}.{assembler}.prodigal.out/{sample}.{assembler}.ffn"),
            uc = os.path.join(
                config["output"]["cobinning"],
                "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.uc.gz")
        output:
            cds = os.path.join(
                config["output"]["cobinning"],
                "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.marker.fa.gz")
        log:
            os.path.join(config["output"]["cobinning"],
                         "logs/marker/{sample}.{assembler}.cds.marker.log")
        run:
            from Bio import SeqIO, bgzf
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

            with open(input.cds, 'r') as ih, bgzf.BgzfWriter(output.cds, 'wb') as oh:
                records = [r for r in SeqIO.parse(ih, "fasta") if r.id in marker_cds]
                cds_count = SeqIO.write(records, oh, "fasta")
                logh.write("Saved %i marker genes to %s\n" % (cds_count, output.cds))
            logh.close()


    rule cobinning_index_cds_marker:
        input:
            os.path.join(
                config["output"]["cobinning"],
                "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.marker.fa.gz")
        output:
            expand("{prefix}.{suffix}",
                   prefix = os.path.join(
                       config["output"]["cobinning"],
                       "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.marker.fa.gz"),
                   suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        log:
            os.path.join(config["output"]["cobinning"],
                         "logs/index/{sample}.{assembler}.cds.marker.index.log")
        params:
            prefix = os.path.join(
                config["output"]["cobinning"],
                "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.marker.fa.gz")
        shell:
            '''
            bowtie2-build \
            {input} \
            {params.prefix} \
            > {log}
            '''


    rule cobinning_get_marker_contigs_coverage:
        input:
            reads = cobinning_input,
            db = expand("{prefix}.{suffix}",
                        prefix = os.path.join(
                            config["output"]["cobinning"],
                            "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.marker.fa.gz"),
                        suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            os.path.join(config["output"]["cobinning"],
                         "coverage/{assembler}/{sample_}/{sample_}.{sample}.{assembler}.coverage.gz")
        log:
            os.path.join(config["output"]["cobinning"],
                         "logs/coverage/{assembler}/{sample_}/{sample_}.{sample}.{assembler}.coverage.log")
        params:
            index = os.path.join(
                config["output"]["cobinning"],
                "cds/{sample}.{assembler}.vsearch.out/{sample}.{assembler}.cds.marker.fa.gz"),
            coverage = os.path.join(
                config["output"]["cobinning"],
                "coverage/{assembler}/{sample_}/{sample_}.{sample}.{assembler}.coverage")
        threads:
            config["params"]["cobinning"]["threads"]
        shell:
            '''
            bowtie2 \
            --threads {threads} \
            -x {params.index} \
            -1 {input.reads[0]} \
            -2 {input.reads[1]} 2> {log} | \
            samtools sort \
            -@{threads} \
            -O BAM - | \
            jgi_summarize_bam_contig_depths \
            --outputDepth {params.coverage} - \
            2>> {log}

            pigz {params.coverage}
            '''


    rule cobinning_all:
        input:
            expand(os.path.join(
                config["output"]["cobinning"],
                "coverage/{assembler}/{sample_}/{sample_}.{sample}.{assembler}.coverage.gz"),
                   assembler=ASSEMBLERS,
                   sample_=SAMPLES.index.unique(),
                   sample=SAMPLES.index.unique()),

            rules.predict_scaftigs_gene_all.input

else:
    rule cobinning_all:
        input:
