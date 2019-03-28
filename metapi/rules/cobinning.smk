rule filter_scaftigs_to_cds:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        scaftigs = expand(os.path.join(config["results"]["cobinning"]["scaftigs"],
                                       "{{sample}}.{{assembler}}.scaftigs.{length}.fa"),
                          length=config["params"]["cobinning"]["scaftigs_length"]),
        cds = os.path.join(config["results"]["cobinning"]["cds"],
                           "{sample}/{sample}.{assembler}.cds.fa"),
        gff = os.path.join(config["results"]["cobinning"]["cds"],
                           "{sample}/{sample}.{assembler}.cds.gff")
    params:
        length = config["params"]["cobinning"]["scaftigs_length"],
        id = "{sample}"
    threads:
        config["params"]["cobinning"]["threads"]
    shell:
        '''
        seqkit seq --min-len {params.length} --threads {threads} --quiet {input} | 
        seqkit replace --pattern "^(.+)" --replacement '{params.id}_$1' |
        tee >(prodigal -d {output.cds} -o {output.gff} -f gff -p meta -q) > {output.scaftigs}
        '''


rule vsearch_clust_cds:
    input:
        cds = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.fa")
    output:
        uc = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.uc")
    params:
        identity = config["params"]["cobinning"]["vsearch"]["identity"]
    log:
        os.path.join(config["logs"]["cobinning"], "{sample}.{assembler}.vsearch_clust.cds.log")
    threads:
        os.path.join(config["cobinning"]["threads"])
    shell:
        '''
        vsearch --threads {threads} --cluster_fast {input.cds} --id {params.identity} --uc {output.uc} 2> {log}
        '''


rule choose_cds_marker:
    input:
        scaftigs = expand(os.path.join(config["results"]["cobinning"]["scaftigs"],
                                       "{{sample}}.{{assembler}}.scaftigs.{length}.fa"),
                          length=config["params"]["cobinning"]["scaftigs_length"]),
        cds = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.fa"),
        uc = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.uc")
    output:
        cds = os.path.join(config["results"]["cobinning"]["cds"], "{sample}/{sample}.{assembler}.cds.marker.fa")
    log:
        os.path.join(config["logs"]["cobinning"], "{sample}.{assembler}.cds.marker.log")
    run:
        from Bio import SeqIO
        import subprocess

        cmd = 'grep -c "^>" %s' % input.scaftigs
        contig_count = int(subprocess.check_output(cmd, shell=True))
        logh = open(log, 'w')

        clust = set()
        hit = set()
        with open(input.rc, 'r') as uch:
            for line in uch:
                item = line.strip().split('\t')
                if 'H' in item[0]:
                    hit.add(item[-1])
                    hit.add(item[-2])
                elif 'C' in item[0]:
                    clust.add(item[-2])

        marker_cds = clust - hit
        marker_contig = set()
        for cds in marker_cds:
            contig_id = '_'.join(cds.split('_')[0:-1])
            marker_contig.add(contig_id)

        marker_ratio = len(marker_contig) / float(contig_count) * 100
        log.write("ratio of contig with marker gene: %.2f%%\n" % marker_ratio)

        with open(input.cds, 'r') as ih, open(output.cds, 'w') as oh:
            records = [r for r in SeqIO.parse(ih, "fasta") if r.id in marker_cds]
            cds_count = SeqIO.write(records, oh, "fasta")
            logh.write("Saved %i marker genes to %s\n" % (cds_count, output.cds))

        logh.close()


rule index_marker_cds:
    input:
        os.path.join(config["results"]["cobinning"]["cds"],
                     "{sample}/{sample}.{assembler}.cds.marker.fa")
    output:
        expand("{prefix}.{suffix}",
               prefix=os.path.join(config["results"]["cobinning"]["cds"],
                                   "{{sample}}/{{sample}}.{{assembler}}.cds.marker.fa"),
               suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    log:
        os.path.join(config["logs"]["cobinning"], "{sample}.{assembler}.cds.marker.index.log")
    params:
        prefix = os.path.join(config["results"]["cobinning"]["cds"],
                              "{sample}/{sample}.{assembler}.cds.marker.fa")
    shell:
        '''
        bowtie2-build {input} {params.prefix} >{log} 2>&1
        '''


'''
rule alignment:
    input:
    output:
    log:
    params:
    shell:
'''





