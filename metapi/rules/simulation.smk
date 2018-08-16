rule genome_download:
    output:
        metadata = os.path.join(config["results"]["simulation"], "species_metadata.tsv")
    params:
        outdir = config["results"]["simulation"],
        taxid = ",".join(config["params"]["simulation"]["taxid"])
    shell:
        '''
        ncbi-genome-download \
            --format fasta,assembly-report \
            --assembly-level complete \
            --taxid {params.taxid} \
            --refseq-category reference \
            --output-folder {params.outdir} \
            --human-readable \
            --retries 3 \
            -m {output.metadata} bacteria
        '''


rule genome_merge:
    input:
        metadata = os.path.join(config["results"]["simulation"], "species_metadata.tsv")
    output:
        expand("{simulation}/{sample}_genome.fa",
               simulation=config["results"]["simulation"],
               sample=_samples.index)
    params:
        outdir = config["results"]["simulation"]
    run:
        import glob
        import gzip
        import re
        import subprocess
        from Bio import SeqIO

        genomes_list = glob.glob(params.outdir + "/refseq/bacteria/*/*.fna.gz")

        def extract_genome(genome_list, out):
            with open(out, 'w') as oh:
                for genome in genome_list:
                    if genome.endswith(".gz"):
                        gh = gzip.open(genome, 'rt')
                    else:
                        gh = open(genome, 'r')
                    for record in SeqIO.parse(gh, 'fasta'):
                        if not re.search(r'plasmid', record.description):
                            SeqIO.write(record, oh, 'fasta')

        extract_genome(genomes_list[0:4], output[0])
        extract_genome(genomes_list[1:5], output[1])
        extract_genome(genomes_list[2:6], output[2])


rule genome_simulate:
    input:
        os.path.join(config["results"]["simulation"], "{sample}_genome.fa")
    output:
        r1 = os.path.join(config["results"]["raw"]["reads"], "{sample}_1.fq.gz"),
        r2 = os.path.join(config["results"]["raw"]["reads"], "{sample}_2.fq.gz"),
        abundance = os.path.join(config["results"]["raw"]["reads"], "{sample}_abundance.txt")
    params:
        model = config["params"]["simulation"]["model"],
        n_reads = config["params"]["simulation"]["n_reads"],
        prefix = os.path.join(config["results"]["raw"]["reads"], "{sample}")
    threads:
        config["params"]["simulation"]["threads"]
    shell:
        '''
        iss generate --cpus {threads} --genomes {input} --n_reads {params.n_reads} --model {params.model} --output {params.prefix}
        pigz {params.prefix}_R1.fastq
        pigz {params.prefix}_R2.fastq
        mv {params.prefix}_R1.fastq.gz {output.r1}
        mv {params.prefix}_R2.fastq.gz {output.r2}
        '''
