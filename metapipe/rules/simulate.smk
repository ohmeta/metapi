# Lactobacillus salivarius UCC118 (firmicutes): 362948
# Lactobacillus paracasei ATCC 334 (firmicutes): 321967
# Klebsiella pneumoniae subsp. pneumoniae HS11286 (enterobacteria): 1125630
# Bifidobacterium bifidum PRL2010 (high GC Gram+): 702459
# Bacteroides thetaiotaomicron VPI-5482 (CFB group bacteria): 226186
# Clostridioides difficile 630 (firmicutes): 272563

rule genome_download:
    input:
        config["results"]["config"]
    output:
        dir = config["results"]["simulate"]["genome"],
        metadata = os.path.join(config["results"]["simulate"]["genome"], "species_metadata.tsv")
    params:
        taxid = ",".join(config["params"]["simulate"]["taxid"])
    shell:
        '''
        echo "downloading Lactobacillus salivarius UCC118 (firmicutes)"
        echo "downloading Lactobacillus paracasei ATCC 334 (firmicutes)"
        echo "downloading Klebsiella pneumoniae subsp. pneumoniae HS11286 (enterobacteria)"
        echo "downloading Bifidobacterium bifidum PRL2010 (high GC Gram+)"
        echo "downloading Bacteroides thetaiotaomicron VPI-5482 (CFB group bacteria)"
        echo "downloading Clostridioides difficile 630 (firmicutes)"
        
        ncbi-genome-download \
            --format fasta,assembly-report \
            --assembly-level complete \
            --taxid {params.taxid} \
            --refseq-category reference \
            --output-folder {output.dir} \
            --human-readable \
            --retries 3 \
            -m {output.metadata} bacteria  
        '''

rule genome_merge:
    input:
        dir = os.path.join(config["results"]["simulate"]["genome"], "refseq/bacteria"),
        metadata = os.path.join(config["results"]["simulate"]["genome"], "species_metadata.tsv")
    output:
        merge_fa = os.path.join(config["results"]["simulate"]["genome"], "merged_genome.fa")
    shell:
        '''
        find {input.dir} -type f -name "*genomic.fna.gz" | xargs -I {} zcat {} >> {output.merge_fa} 
        echo "hello"
        '''

'''
|>NC_014638.1 Bifidobacterium bifidum PRL2010 chromosome, complete genome                                                      
│>NC_004663.1 Bacteroides thetaiotaomicron VPI-5482 chromosome, complete genome                                                
│>NC_004703.1 Bacteroides thetaiotaomicron VPI-5482 plasmid p5482, complete sequence                                           
│>NC_008526.1 Lactobacillus paracasei ATCC 334 chromosome, complete genome                                                     
│>NC_008502.1 Lactobacillus paracasei ATCC 334 plasmid 1, complete sequence                                                    
│>NC_009089.1 Clostridioides difficile 630, complete genome
│>NC_008226.2 Clostrididioides difficile 630 plasmid pCD630, complete replicon                                                 
│>NC_007929.1 Lactobacillus salivarius UCC118 chromosome, complete genome                                                      
│>NC_007930.1 Lactobacillus salivarius UCC118 plasmid pMP118, complete sequence                                                
│>NC_006529.1 Lactobacillus salivarius UCC118 plasmid pSF118-20, complete sequence                                             
│>NC_006530.1 Lactobacillus salivarius UCC118 plasmid pSF118-44, complete sequence                                             
│>NC_016845.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 chromosome, complete genome                                      
│>NC_016838.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 plasmid pKPHS1, complete sequence                                
│>NC_016846.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 plasmid pKPHS2, complete sequence                                
│>NC_016839.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 plasmid pKPHS3, complete sequence                                
│>NC_016840.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 plasmid pKPHS4, complete sequence                                
│>NC_016847.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 plasmid pKPHS5, complete sequence                                
│>NC_016841.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 plasmid pKPHS6, complete sequence
'''

rule simulate:
    input:
        os.path.join(config["results"]["simulate"]["genome"], "merged_genome.fa")
    output:
        r1 = os.path.join(config["results"]["simulate"]["genome"], "simulate_100X")
    params:
        model = config["params"]["simulate"]["model"],
        n_genomes =  len(config["params"]["simulate"]["taxid"]),
        n_reads = config["params"]["simulate"]["coverage"]["X100"],
        output_prefix = config["params"]["simulate"]["output_prefix"]["X100"]
    shell:
        '''
        iss generate 
            --cpus 8 \
            --genomes {intput} \
            --n_genomes {params.n_genomes} \
            --n_reads {params.n_reads} \
            --model {params.model} \
            --output {params.output_prefix}
        '''