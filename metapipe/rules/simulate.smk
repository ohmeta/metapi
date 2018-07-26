rule genome_download:
    output:
        metadata = os.path.join(config["results"]["simulate"]["genome"], "species_metadata.tsv")
    params:
        dir = config["results"]["simulate"]["genome"],
        taxid = ",".join(config["params"]["simulate"]["taxid"])
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
        metadata = os.path.join(config["results"]["simulate"]["genome"], "species_metadata.tsv")
    output:
        merge_fa = os.path.join(config["results"]["simulate"]["genome"], "merged_genome.fasta")
    params:
        dir = os.path.join(config["results"]["simulate"]["genome"])
    shell:
        '''
        zcat {params.dir}/refseq/bacteria/*/*genomic.fna.gz > {output.merge_fa}
        '''

rule simulate:
    input:
        os.path.join(config["results"]["simulate"]["genome"], "merged_genome.fasta")
    output:
        r1 = os.path.join(config["results"]["simulate"]["genome"], "simulate_100X_R1.fastq"),
        r2 = os.path.join(config["results"]["simulate"]["genome"], "simulate_100X_R2.fastq"),
        abundance = os.path.join(config["results"]["simulate"]["genome"], "simulate_100X_abundance.txt")
    params:
        model = config["params"]["simulate"]["model"],
        n_genomes =  len(config["params"]["simulate"]["taxid"]),
        n_reads = config["params"]["simulate"]["coverage"]["X100"],
        output_prefix = os.path.join(config["results"]["simulate"]["genome"],
                                     config["params"]["simulate"]["output_prefix"]["X100"])
    threads:
        config["params"]["simulate"]["threads"]
    shell:
        '''
        iss generate --cpus {threads} --genomes {input} --n_genomes {params.n_genomes} --n_reads {params.n_reads} --model {params.model} --output {params.output_prefix}
        '''
