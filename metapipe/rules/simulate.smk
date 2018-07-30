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
        r1 = os.path.join(config["results"]["simulation"]["genomes"], "simulate_100X.1.fastq"),
        r2 = os.path.join(config["results"]["simulation"]["genomes"], config["results"]["simulation"]["out_prefix"]".2.fastq"),
        abundance = os.path.join(config["results"]["simulation"]["genomes"], "simulate_100X_abundance.txt")
    params:
        model = config["params"]["simulation"]["model"],
        n_genomes = len(config["params"]["simulation"]["taxid"]),
        n_reads = config["params"]["simulation"]["n_reads"],
        output_prefix = os.path.join(config["results"]["simulation"]["genomes"],
                                     config["params"]["simulation"]["output_prefix"])
    threads:
        config["params"]["simulate"]["threads"]
    shell:
        '''
        iss generate --cpus {threads} --genomes {input} --n_genomes {params.n_genomes} --n_reads {params.n_reads} --model {params.model} --output {params.output_prefix}
        mv {config["results"]["simulate"]["genome"]/simulate_100X_R1.fastq} {output.r1}
        mv {config["results"]["simulate"]["genome"]/simulate_100X_R2.fastq} {output.r2}
        '''
