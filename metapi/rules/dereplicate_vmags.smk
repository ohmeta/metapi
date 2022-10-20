# reference 1: MGV
# https://github.com/snayfach/MGV/tree/master/ani_cluster

# reference 2: CheckV
# https://bitbucket.org/berkeleylab/checkv/src/master/
# Supporting code

# Rapid genome clustering based on pairwise ANI

# First, create a blast+ database:
# makeblastdb -in <my_seqs.fna> -dbtype nucl -out <my_db>

# Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
# blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -o <my_blast.tsv> -num_threads 32

# Next, calculate pairwise ANI by combining local alignments between sequence pairs:
# anicalc.py -i <my_blast.tsv> -o <my_ani.tsv>

# Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
# aniclust.py --fna <my_seqs.fna> --ani <my_ani.tsv> --out <my_clusters.tsv> --min_ani 95 --min_tcov 85 --min_qcov 0 


CLUSTER_VMAGS_GROUPS = CHECKV_GROUPS.loc[:, ["binning_group", "assembly_group", "identifier"]].drop_duplicates()


rule dereplicate_vmags_prepare:
    input:
        expand(os.path.join(config["output"]["check"],
               "data/checkv/{binning_group}.{assembly_group}.{{assembler}}/{identifier}/vMAG_hmq.fa"),
                zip,
                binning_group=CLUSTER_VMAGS_GROUPS["binning_group"],
                assembly_group=CLUSTER_VMAGS_GROUPS["assembly_group"],
                identifier=CLUSTER_VMAGS_GROUPS["identifier"])
    output:
        os.path.join(config["output"]["dereplicate"], "genomes/virome/input/vMAGs_hmq.{assembler}.fa")
    shell:
        '''
        cat {input} > {output} 
        '''


rule dereplicate_vmags_build_db:
    input:
        os.path.join(config["output"]["dereplicate"], "genomes/virome/input/vMAGs_hmq.{assembler}.fa")
    output:
        os.path.join(config["output"]["dereplicate"], "genomes/virome/blastdb/vMAGs_hmq.{assembler}.blastdb")
    benchmark:
        os.path.join(config["output"]["dereplicate"], "benchmark/blastdb/dereplicate_vmags_build_db.{assembler}.benchmark.txt")
    log:
        os.path.join(config["output"]["dereplicate"], "logs/blastdb/dereplicate_vmags_build_db.{assembler}.blastdb.log")
    conda:
        config["envs"]["blast"]
    threads:
        config["params"]["dereplicate_vmags"]["threads"]
    shell:
        '''
        makeblastdb \
        -in {input} \
        -dbtype nucl \
        -out {output} \
        >{log} 2>&1
        '''


rule dereplicate_vmags_blastn:
    input:
        fa = os.path.join(config["output"]["dereplicate"],
                          "genomes/virome/input/vMAGs_hmq.{assembler}.fa"),
        db = os.path.join(config["output"]["dereplicate"],
                          "genomes/virome/blastdb/vMAGs_hmq.{assembler}.blastdb")
    output:
        os.path.join(config["output"]["dereplicate"],
                     "genomes/virome/blastout/vMAGs_hmq.{assembler}.blast.fmt6.tsv")
    benchmark:
        os.path.join(config["output"]["dereplicate"],
                     "benchmark/blast/dereplicate_vmags_blastn.{assembler}.benchmark.txt")
    log:
        os.path.join(config["output"]["dereplicate"],
                     "logs/blast/dereplicate_vmags_blastn.{assembler}.log")
    conda:
        config["envs"]["blast"]
    params:
        max_target_seqs = config["params"]["dereplicate_vmags"]["blast"]["max_target_seqs"],
        prec_identity = config["params"]["dereplicate_vmags"]["blast"]["prec_identity"]
    threads:
        config["params"]["dereplicate_vmags"]["threads"]
    shell:
        '''
        blastn \
        -num_threads {threads} \
        -query {input.fa} \
        -db {input.db} \
        -out {output} \
        -outfmt '6 std qlen slen' \
        -max_target_seqs {params.max_target_seqs} \
        -perc_identity {params.prec_identity} \
        >{log} 2>&1
        '''
 

rule dereplicate_vmags_compute_ani:
    input:
        os.path.join(config["output"]["dereplicate"],
                     "genomes/virome/blastout/vMAGs_hmq.{assembler}.blast.fmt6.tsv")
    output:
        os.path.join(config["output"]["dereplicate"],
                     "genomes/virome/ani/vMAGs_hmq.{assembler}.ani.tsv")
    benchmark:
        os.path.join(config["output"]["dereplicate"],
                     "benchmark/virome/dereplicate_vmags_compute_ani.{assembler}.benchmark.txt")
    log:
        os.path.join(config["output"]["dereplicate"],
                     "logs/virome/dereplicate_vmags_compute_ani.{assembler}.log")
    params:
        anicalc = os.path.join(WRAPPER_DIR, "anicalc.py")
    threads:
        1
    shell:
        '''
        python {params.anicalc} \
        -i {input} -o {output} \
        >{log} 2>&1
        '''


rule dereplicate_vmags_clust:
    input:
        fa = os.path.join(config["output"]["dereplicate"],
                          "genomes/virome/input/vMAGs_hmq.{assembler}.fa"),
        ani = os.path.join(config["output"]["dereplicate"],
                           "genomes/virome/ani/vMAGs_hmq.{assembler}.ani.tsv")
    output:
        os.path.join(config["output"]["dereplicate"],
                     "genomes/virome/cluster/vMAGs_hmq.{assembler}.cluster.tsv")
    benchmark:
        os.path.join(config["output"]["dereplicate"],
                     "benchmark/virome/dereplicate_vmags_clust.{assembler}.benchmark.txt")
    log:
        os.path.join(config["output"]["dereplicate"],
                     "logs/virome/dereplicate_vmags_clust.{assembler}.log")
    params:
        aniclust = os.path.join(WRAPPER_DIR, "aniclust.py"),
        min_ani = config["params"]["dereplicate_vmags"]["blast"]["min_ani"],
        min_tcov = config["params"]["dereplicate_vmags"]["blast"]["min_tcov"],
        min_qcov = config["params"]["dereplicate_vmags"]["blast"]["min_qcov"]
    threads:
        1
    shell:
        '''
        python {params.aniclust} \
        --fna {input.fa} \
        --ani {input.ani} \
        --out <output> \
        --min_ani 95 \
        --min_tcov 85 \
        --min_qcov 0 \
        >{log} 2>&1
        '''


#rule dereplicate_vmags:
#    input:
#    output:
#    shell:


rule dereplicate_vmags_all:
    input:
        expand(
            os.path.join(
                config["output"]["dereplicate"],
                "genomes/virome/cluster/vMAGs_hmq.{assembler}.cluster.tsv"),
                assembler=ASSEMBLERS)
 

rule dereplicate_all:
    input:
        rules.dereplicate_gene_all.input,
        rules.dereplicate_mags_all.input,
        rules.dereplicate_vmags_all.input


localrules:
    dereplicate_vmags_all
    dereplicate_all