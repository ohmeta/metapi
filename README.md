# **metapipe**  (work in process)

hello, metagenomics!

## brother project

* [sunbeam](https://github.com/sunbeam-labs/sunbeam)
* [atlas](https://github.com/pnnl/atlas)

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
  * [workflow language definition](https://software.broadinstitute.org/wdl/)
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
        include: "rules/trim.smk"
        include: "rules/rmhost.smk"
        include: "rules/qcreport.smk"
        include: "rules/assembly.smk"
        include: "rules/alignment.smk"
        include: "rules/binning.smk"
        include: "rules/checkm.smk"
        include: "rules/drep.smk"
        include: "rules/annotation.smk"
        include: "rules/reassembly.smk"
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

* install dependencies
  * [pigz](https://zlib.net/pigz/)
  * [sickle](https://github.com/najoshi/sickle)
  * [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
  * [bwa](https://github.com/lh3/bwa)
  * [samtools](https://github.com/samtools/samtools)
  * [spades](https://github.com/ablab/spades)
  * [idba](https://github.com/loneknightpy/idba)
  * [megahit](https://github.com/voutcn/megahit)
  * [metabat](https://bitbucket.org/berkeleylab/metabat)
  * [checkm](https://github.com/Ecogenomics/CheckM)
  * [mash](https://github.com/marbl/Mash)
  * [sourmash](https://github.com/dib-lab/sourmash)
  * [centrifuge](https://github.com/infphilo/centrifuge)
  * [prokka](https://github.com/tseemann/prokka)

  ```bash
  # in python3 environment
  conda install pigz snakemake sickle-trim bbmap bwa samtools spades idba megahit mash sourmash centrifuge prokka
  conda install -c ursky metabat2
  
  # in python2 envrionment
  conda install checkm-genome
  
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
    git clone https://github.com/ohmeta/metapipe
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
    params:
        bwa_t = 8,
        samtools_t = 8
    shell:
        "bwa mem -t {params.bwa_t} {input.prefix} {input.r1} {input.r2} | "
        "samtools view -@{params.samtools_t} -hbS - | "
        "tee >(samtools flagstat -@{params.samtools_t} - > {output.flagstat}) | "
        "samtools sort -@{params.samtools_t} -o {output.bam} -"
    ```

* a simulated metagenomics data test(uncomplete)

    ```bash
    # in metapipe/.test directory
    cd .test

    # look
    snakemake --dag | dot -Tsvg > dat.svg
    ```
    <img src=".test/dat.svg">

    ```bash
    # run
    snakemake
    ```

* a real world metagenomics data process(uncomplete)

    ```bash
    # in metapipe directory
    # look
    snakemake --dag | dot -Tsvg > dat.svg
    ```
    <img src="dat.svg">

    ```bash
    # run
    snakemake
    ```

## reference

## contributor

* [alienzj](https://github.com/alienzj)
* [Margaret Chu](https://github.com/magcurly)