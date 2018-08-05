rule genome_download:
    output:
        metadata = os.path.join(config["results"]["simulation"]["genomes"], "species_metadata.tsv")
    params:
        dir = config["results"]["simulation"]["genomes"],
        taxid = ",".join(config["params"]["simulation"]["taxid"])
    shell:
        '''
        ncbi-genome-download \
            --format fasta,assembly-report \
            --assembly-level complete \
            --taxid {params.taxid} \
            --refseq-category reference \
            --output-folder {params.dir} \
            --human-readable \
            --retries 3 \
            -m {output.metadata} bacteria
        '''

rule genome_merge:
    input:
        metadata = os.path.join(config["results"]["simulation"]["genomes"], "species_metadata.tsv")
    output:
        merge_fa = os.path.join(config["results"]["simulation"]["genomes"], "merged_genome.fasta")
    params:
        dir = os.path.join(config["results"]["simulation"]["genomes"])
    shell:
        '''
        zcat {params.dir}/refseq/bacteria/*/*genomic.fna.gz > {output.merge_fa}
        '''

rule genome_simulate:
    input:
        os.path.join(config["results"]["simulation"]["genomes"], "merged_genome.fasta")
    output:
        r1 = os.path.join(config["results"]["simulation"]["genomes"], config["params"]["simulation"]["output_prefix"] + "_1.fq.gz"),
        r2 = os.path.join(config["results"]["simulation"]["genomes"], config["params"]["simulation"]["output_prefix"] + "_2.fq.gz"),
        abundance = os.path.join(config["results"]["simulation"]["genomes"], config["params"]["simulation"]["output_prefix"] + "_abundance.txt")
    params:
        model = config["params"]["simulation"]["model"],
        n_genomes = config["params"]["simulation"]["n_genomes"],
        n_reads = config["params"]["simulation"]["n_reads"],
        prefix = os.path.join(config["results"]["simulation"]["genomes"],
                                     config["params"]["simulation"]["output_prefix"])
    threads:
        config["params"]["simulation"]["threads"]
    shell:
        '''
        iss generate --cpus {threads} --genomes {input} --n_genomes {params.n_genomes} --n_reads {params.n_reads} --model {params.model} --output {params.prefix}
        pigz {params.prefix}_R1.fastq
        pigz {params.prefix}_R2.fastq
        mv {params.prefix}_R1.fastq.gz {output.r1}
        mv {params.prefix}_R2.fastq.gz {output.r2}
        '''
