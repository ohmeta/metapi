# **metapi**

hello, metagenomics!

## brother project

* [atlas](https://github.com/pnnl/atlas)
* [sunbeam](https://github.com/sunbeam-labs/sunbeam)

## motivation

  we all need a metagenomics pipeline for academic research.

## principle

* bind intelligense together
  * [github](https://github.com/search?q=metagenomics)
  * why we here?
* do not make wheels
  * [biopython](https://github.com/biopython/biopython)
  * [bioperl](http://bioperl.org)
  * [seqan](https://github.com/seqan/seqan)
* make full use of pipeline execution engine
  * [make](https://www.gnu.org/software/make/manual/make.html)
  * [snakemake](https://bitbucket.org/snakemake/snakemake)
  * [common workflow language](https://github.com/common-workflow-language/common-workflow-language)
  * [workflow definition language](https://software.broadinstitute.org/wdl/)
* make full use of awesome bioinformatics tools
  * [bwa](https://github.com/lh3/bwa)
  * [samtools](https://github.com/samtools/samtools)
  * [spades](https://github.com/ablab/spades)
  * [idba](https://github.com/loneknightpy/idba)
  * [megahit](https://github.com/voutcn/megahit)
  * [metabat](https://bitbucket.org/berkeleylab/metabat)
  * [checkm](https://github.com/Ecogenomics/CheckM)
* robust and module, extensible, update
  * one rule, one module
  * one module, one analysis
* welcome to PR

## design

* execution module
    ```python
    # Snakefile
        include: "rules/step.smk"
        include: "rules/simulation.smk"
        include: "rules/fastqc.smk"
        include: "rules/trimming.smk"
        include: "rules/rmhost.smk"
        include: "rules/assembly.smk"
        include: "rules/alignment.smk"
        include: "rules/binning.smk"
        include: "rules/checkm.smk"
        include: "rules/dereplication.smk"
        include: "rules/classification.smk"
        include: "rules/annotation.smk"
        include: "rules/profilling.smk"
    ```

* analysis module
  * raw data report
  * quality control
  * remove host sequences
  * assembly
  * assembly evaluation
  * binning
  * checkm
  * dereplication
  * bins profile
  * taxonomy classification
  * genome annotation
  * function annotation

* test module
  * execution test
  * analysis test

## install

* install dependencies*
  * [snakemake](https://snakemake.readthedocs.io)
  * [pigz](https://zlib.net/pigz/)
  * [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
  * [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)
  * [OAFilter](https://github.com/Scelta/OAFilter)
  * [sickle](https://github.com/najoshi/sickle)
  * [fastp](https://github.com/OpenGene/fastp)
  * [MultiQC](https://github.com/ewels/MultiQC)
  * [bwa](https://github.com/lh3/bwa)
  * [samtools](https://github.com/samtools/samtools)
  * [bbmap](https://sourceforge.net/projects/bbmap)
  * [spades](https://github.com/ablab/spades)
  * [idba](https://github.com/loneknightpy/idba)
  * [megahit](https://github.com/voutcn/megahit)
  * [quast](https://sourceforge.net/projects/quast/)
  * [MetaBat](https://bitbucket.org/berkeleylab/metabat)
  * [MaxBin2](http://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html)
  * [CheckM](https://github.com/Ecogenomics/CheckM)
  * [drep](https://github.com/MrOlm/drep)
  * [prokka](https://github.com/tseemann/prokka)
  * [metaphlan2](https://bitbucket.org/biobakery/metaphlan2)

  ```bash
  # in python3 environment
  conda install snakemake pigz ncbi-genome-download sickle-trim fastp bwa samtools bbmap spades idba megahit maxbin2 prokka
  conda install -c ursky metabat2
  pip install drep insilicoseq

  # in python2 envrionment
  conda install quast checkm-genome metaphlan2

  # database configuration
  wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
  mkdir checkm_data
  cd checkm_data
  tar -xzvf ../checkm_data_2015_01_16.tar.gz
  cd ..
  ln -s checkm_data checkm_data_latest

  # activate python2 environment where checkm in
  checkm data setRoot checkm_data_latest
  ```

* install metapipe

    ```bash
    git clone https://github.com/ohmeta/metapi
    # or
    pip install metapi
    ```

## example

* snakemake了解一下:)

    ```python
    rule bwa_mem:
    input:
        r1 = "fastq/sample_1.fq.gz",
        r2 = "fastq/sample_2.fq.gz",
        ref = "ref/ref.index
    output:
        bam = "sample.sort.bam",
        stat = "sample_flagstat.txt"
    threads:
        8
    shell:
        '''
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools view -@{threads} -hbS - | \
        tee >(samtools flagstat -@{threads} - > {output.stat}) | \
        samtools sort -@{threads} -o {output.bam} -
        '''
    ```

* a simulated metagenomics data test(uncomplete)

    ```bash
    # in metapi/example/basic_test directory
    cd example/basic_test

    # look
    snakemake --dag | dot -Tsvg > dat.svg
    ```
    <img src="examples/basic_test/dat.svg">

    ```bash
    # run on local
    snakemake

    # run on SGE cluster
    snakemake \
    --jobs 80 \
    --cluster "qsub -S /bin/bash -cwd -q {queue} -P {project} -l vf={mem},p={cores} -binding linear:{cores}"
    ```

* a real world metagenomics data process(uncomplete)

    ```bash
    # in metapipe directory
    # look
    cd metapi
    snakemake --dag | dot -Tsvg > ../docs/dat.svg
    ```
    <img src="docs/dat.svg">

    ```bash
    # run on local
    snakemake \
    --cores 8 \
    --snakefile metapi/Snakefile \
    --configfile metapi/metaconfig.yaml \
    --until all

    # run on SGE cluster
    snakemake \
    --snakefile metapi/Snakefile \
    --configfile metapi/metaconfig.yaml \
    --cluster-config metapi/metacluster.yaml \
    --jobs 80 \
    --cluster "qsub -S /bin/bash -cwd -q {cluster.queue} -P {cluster.project} -l vf={cluster.mem},p={cluster.cores} -binding linear:{cluster.cores} -o {cluster.output} -e {cluster.error}"
    --latency-wait 360 \
    --until all
    ```
    
    ```bash
    metapi --help
    
    usage: metapi [subcommand] [options]
    
    metapi, a metagenomics data process pipeline

    optional arguments:
        -h, --help     show this help message and exit
        -v, --version  print software version and exit

    available subcommands:
    
        init         a metagenomics project initialization
        simulation   a simulation on metagenomics data
        workflow     a workflow on real metagenomics data

    ```
    
    ```bash
    metapi init --help
    
    usage: metapi init [-h] [-d <str>] [-q QUEUE] [-p PROJECT] [-s SAMPLES]
                       [-b {raw,assembly}] [-a [<str> [<str> ...]]]

    a metagenomics project initialization
    
    arguments:
        -h, --help            show this help message and exit
        -d <str>, --workdir <str>
                              project workdir
        -q QUEUE, --queue QUEUE
                          　　　cluster queue
        -p PROJECT, --project PROJECT
                          　　　project id
    　　　-s SAMPLES, --samples SAMPLES
                          　　　raw fastq samples list
    　　　-b {raw,assembly}, --begin {raw,assembly}
                          　　　begin to run pipeline from a specific step
    　　　-a [<str> [<str> ...]], --assembler [<str> [<str> ...]]
                          　　　support metaspades, idba_ud, megahit
    ```
    
    ```bash
    metapi simulation --help
    
    usage: metapi simulation [-h] [-d <str>] [-t [<int> [<int> ...]]]
                             [-g <genomes.fasta>] [-ng <int>] [-nr <str>]
                             [-m {hiseq,novaseq,miseq}]
                             [-u {genome_download,
                                  genome_merge,
                                  genome_simulate,
                                  fastqc,multiqc_fastqc,
                                  trimming_fastp,
                                  multiqc_fastp,
                                  build_host_index,
                                  rmhost,
                                  assembly_megahit,
                                  assembly_idba_ud,
                                  assembly_metaspades,
                                  coassembly_megahit,
                                  metaquast_megahit,
                                  multiqc_metaquast,
                                  build_asmfa_index,
                                  align_reads_to_asmfa,
                                  coverage_metabat2,
                                  coverage_maxbin2,
                                  binning_metabat2,
                                  binning_maxbin2,
                                  checkm_lineage_wf,
                                  prokka_bins,
                                  metaphlan2_profilling,
                                  metaphlan2_merge,
                                  burst_reads,all}]
　　
   a simulation on metagenomics data
    
   arguments:
        -h, --help            show this help message and exit
        -d <str>, --workdir <str>
                              project workdir
        -t [<int> [<int> ...]], --taxid [<int> [<int> ...]]
                            　reference database species id(sapce-separated)
　　　　　-g <genomes.fasta>, --genomes <genomes.fasta>
　　　　　　　　　　　　　　　　　　genomes fasta, default: None
　　　　　-ng <int>, --n_genomes <int>
　　　　　　　　　　　　　　　　　　genomes number, default: 6
        -nr <str>, --n_reads <str>
                              reads coverage, default: 5M
        -m {hiseq,novaseq,miseq}, --model {hiseq,novaseq,miseq}
                              reads error model, default: hiseq
        -u {genome_download,
            genome_merge,
            genome_simulate,
            fastqc,
            multiqc_fastqc,
            trimming_fastp,
            multiqc_fastp,
            build_host_index,
            rmhost,
            assembly_megahit,
            assembly_idba_ud,
            assembly_metaspades,
            coassembly_megahit,
            metaquast_megahit,
            multiqc_metaquast,
            build_asmfa_index,
            align_reads_to_asmfa,
            coverage_metabat2,
            coverage_maxbin2,
            binning_metabat2,
            binning_maxbin2,
            checkm_lineage_wf,
            prokka_bins,
            metaphlan2_profilling,
            metaphlan2_merge,
            burst_reads,all}, 
        --step {genome_download,
                genome_merge,
                genome_simulate,
                fastqc,
                multiqc_fastqc,
                trimming_fastp,
                multiqc_fastp,
                build_host_index,
                rmhost,
                assembly_megahit,
                assembly_idba_ud,
                assembly_metaspades,
                coassembly_megahit,
                metaquast_megahit,
                multiqc_metaquast,
                build_asmfa_index,
                align_reads_to_asmfa,
                coverage_metabat2,
                coverage_maxbin2,
                binning_metabat2,
                binning_maxbin2,
                checkm_lineage_wf,
                prokka_bins,
                metaphlan2_profilling,
                metaphlan2_merge,
                burst_reads,all}
                run step
    ```
