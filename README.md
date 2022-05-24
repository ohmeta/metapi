[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![ohmeta-badge](https://img.shields.io/badge/install%20with-ohmeta-brightgreen.svg?style=flat)](http://anaconda.org/ohmeta)
[![PyPI version](https://badge.fury.io/py/metapi.svg)](https://badge.fury.io/py/metapi)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metapi/badges/downloads.svg)](https://anaconda.org/bioconda/metapi)

<div align=center><img width="500" height="280" src="docs/logo.svg"/></div>

A general metagenomics data mining system focus on robust microbiome research.

## Overview
### MAG workflow (WIP figure)
<div align=center><img width="600" height="800" src="docs/mag_workflow.svg"></div>

## Installation

metapi works with Python 3.7+.
You can install it via [bioconda](https://bioconda.github.io/):

```
➤ mamba install -c conda-forge -c bioconda metapi

# It is recommended to install the latest version
➤ mamba install -c conda-forge -c bioconda metapi=2.1.1
```

Or via pip:

```
➤ pip3 install metapi=2.1.1
```

## Run

### Help

```
➤ metapi --help

  .___  ___.  _______ .___________.    ___      .______    __
  |   \/   | |   ____||           |   /   \     |   _  \  |  |
  |  \  /  | |  |__   `---|  |----`  /  ^  \    |  |_)  | |  |
  |  |\/|  | |   __|      |  |      /  /_\  \   |   ___/  |  |
  |  |  |  | |  |____     |  |     /  _____  \  |  |      |  |
  |__|  |__| |_______|    |__|    /__/     \__\ | _|      |__|

            Omics for All, Open Source for All

 A general metagenomics data mining system focus on robust microbiome research.

    optional arguments:
    -h, --help     show this help message and exit
    -v, --version  print software version and exit

    available subcommands:

    init         init project
    mag_wf       metagenome-assembly-genome pipeline
    gene_wf      metagenome-assembly-gene pipeline
    sync         metapi sync project
```

### Init

```
➤ metapi init --help
usage: metapi init [-h] [-d WORKDIR] [--check-samples] [-s SAMPLES]
                        [-b {simulate,trimming,rmhost,assembly,binning,checkm}]
                        [--trimmer {oas1,sickle,fastp}]
                        [--rmhoster {soap,bwa,bowtie2,minimap2,kraken2,kneaddata}]
                        [--assembler {idba-ud,megahit,metaspades,spades,opera-ms} [{idba-ud,megahit,metaspades,spades,opera-ms} ...]]
                        [--binner BINNER [BINNER ...]]

optional arguments:
  -h, --help            show this help message and exit
  -d, --workdir WORKDIR
                        project workdir (default: ./)
  --check-samples       check samples, default: False
  -s, --samples SAMPLES
                        desired input:
                        samples list, tsv format required.
                        
                        if begin from trimming, rmhost, or assembly:
                            if it is fastq:
                                the header is: [sample_id, assembly_group, binning_group, fq1, fq2]
                            if it is sra:
                                the header is: [sample_id, assembly_group, binning_group, sra]
                        
                        if begin from simulate:
                                the header is: [id, genome, abundance, reads_num, model]
                        
  -b, --begin {simulate,trimming,rmhost,assembly,binning,checkm}
                        pipeline starting point (default: trimming)
  --trimmer {oas1,sickle,fastp}
                        which trimmer used (default: fastp)
  --rmhoster {soap,bwa,bowtie2,minimap2,kraken2,kneaddata}
                        which rmhoster used (default: bowtie2)
  --assembler {idba-ud,megahit,metaspades,spades,opera-ms} [{idba-ud,megahit,metaspades,spades,opera-ms} ...]
                        which assembler used, required when begin with binning, can be changed in config.yaml (default: ['metaspades'])
  --binner BINNER [BINNER ...]
                        wchich binner used (default: ['metabat2', 'concoct', 'maxbin2', 'vamb', 'dastools'])
```
- **Note**  
  * When we do `metapi init`, metapi will help us to create project structure,
  include `config.yaml`, `profiles/`, `envs/`, `logs/` and `results/`.
    - `config.yaml`: workflow configuration, can be edited
    - `profiles/`: used when run pipeline on cluster, can be edited
      + For `SGE` user, please see `profiles/sge/cluster.yaml`
      + For `SLURM` user, please see `profiles/slurm/cluster.yaml`
    - `envs/`: conda environment file for some rules, can be edited to update software version by yourself

  * We recommand you to create a directory, and `cd` it and put `samples.tsv` in the directory, then do ```metapi init -s samples.tsv -d ./```.

  * Samples file is required when init project. Samples input format requirement can see [here](#inputtag).

  * Begin point can be `simulate`, `trimming`, `rmhost`, `assembly` and `binning`, it means:
    - simulate: simulate reads
    - trimming: samples need trimming, do quality control
    - rmhost: samples is quality control well, only need to remove host sequence, will not do trimming
    - assembly: samples is clean, just do assembly, will not do trimming and rmhost
    -  (WIP) binning: supply samples and assembly, then do binning, will no do trimming, rmhost and assembly

  * When metapi init executed, then corresponding configuration will be writen into `config.yaml`.  
    Of course, you can edit `config.yaml` to update config, then when run `metapi mag_wf` or `metapi gene_wf`,  
    `metapi` will understand it. Just edit it, see what will happen.


### mag_wf

```
➤ metapi mag_wf --help
          [--cluster-engine {slurm,sge,lsf,pbs-torque}]
          [--wait WAIT] [--use-conda] [--conda-prefix CONDA_PREFIX]
          [--conda-create-envs-only]
          [TASK]

positional arguments:
  TASK                  pipeline end point. Allowed values are
  
  simulate_all,
  prepare_short_reads_all,
  prepare_long_reads_all,
  prepare_reads_all,
  raw_fastqc_all,
  raw_report_all, raw_all,
  trimming_oas1_all,
  trimming_sickle_all,
  trimming_fastp_all,
  trimming_report_all,
  trimming_all,
  rmhost_soap_all,
  rmhost_bwa_all,
  rmhost_bowtie2_all,
  rmhost_minimap2_all,
  rmhost_kraken2_all,
  rmhost_kneaddata_all,
  rmhost_report_all,
  rmhost_all,
  qcreport_all,
  assembly_megahit_all,
  assembly_idba_ud_all,
  assembly_metaspades_all,
  assembly_spades_all,
  assembly_plass_all,
  assembly_opera_ms_all, 
  assembly_metaquast_all,
  assembly_report_all,
  assembly_all,
  alignment_base_depth_all,
  alignment_report_all,
  alignment_all,
  binning_metabat2_coverage_all,
  binning_metabat2_all,
  binning_maxbin2_all, 
  binning_concoct_all,
  binning_graphbin2_all,
  binning_dastools_all,
  binning_vamb_prepare_all, 
  binning_vamb_all,
  binning_report_all,
  binning_all,
  identify_virsorter2_all,
  identify_deepvirfinder_all,
  identify_single_all,
  identify_phamb_all,
  identify_multi_all,
  identify_all,
  predict_scaftigs_gene_prodigal_all,
  predict_scaftigs_gene_prokka_all,
  predict_bins_gene_prodigal_all,
  predict_bins_gene_prokka_all,
  predict_scaftigs_gene_all,
  predict_bins_gene_all,
  predict_all,
  checkm_all,
  dereplicate_mags_drep_all,
  dereplicate_mags_all,
  taxonomic_all,
  upload_sequencing_all,
  upload_assembly_all,
  upload_all,
  all (default: all)

optional arguments:
  -h, --help            show this help message and exit
  -d, --workdir WORKDIR
                        project workdir (default: ./)
  --check-samples       check samples, default: False
  --config CONFIG       config.yaml (default: ./config.yaml)
  --profile PROFILE     cluster profile name (default: ./profiles/slurm)
  --cores CORES         all job cores, available on '--run-local' (default: 240)
  --local-cores LOCAL_CORES
                        local job cores, available on '--run-remote' (default: 8)
  --jobs JOBS           cluster job numbers, available on '--run-remote' (default: 30)
  --list                list pipeline rules
  --debug               debug pipeline
  --dry-run             dry run pipeline
  --run-local           run pipeline on local computer
  --run-remote          run pipeline on remote cluster
  --cluster-engine {slurm,sge,lsf,pbs-torque}
                        cluster workflow manager engine, support slurm(sbatch) and sge(qsub) (default: slurm)
  --wait WAIT           wait given seconds (default: 60)
  --use-conda           use conda environment
  --conda-prefix CONDA_PREFIX
                        conda environment prefix (default: ~/.conda/envs)
  --conda-create-envs-only
                        conda create environments only
```

### mag_wf example

```
# init project
➤ metapi init -d . -s samples.tsv -b trimming

# create conda environments, which will create envs at ~/.conda/envs
➤ metapi mag_wf --conda-create-envs-only

# run pipeline with conda
➤ metapi mag_wf all --use-conda --run-local

# run raw_fastqc
➤ metapi mag_wf raw_fastqc_all --run-local

# run trimming
➤ metapi mag_wf trimming_all --run-local

# run rmhost
➤ metapi mag_wf rmhost_all --run-local

# run qc report
➤ metapi mag_wf qcreport_all --run-local

# run assembly
➤ metapi mag_wf assembly_all --run-local

# run binning (microbial MAG)
➤ metapi mag_wf binning_all --run-local

# run identify (virus MAG)
➤ metapi mag_wf identify_all --run-local

# run gene predict
➤ metapi mag_wf predict_all --run-local

# run checkm
➤ metapi mag_wf checkm_all --run-local

#run taxonomic 
➤ metapi mag_wf taxonomic_all --run-local

# run mag_wf all
➤ metapi mag_wf --run-local

# run mag_wf all on SGE/SLURM cluster
➤ metapi mag_wf --run-remote

# run mag_wf all on SGE/SLURM cluster using conda
➤ metapi mag_wf --run-remote --use-conda
```

## Input details
<a name="inputtag"></a>

+ The input fastq need gzip compress.

+ If `sample_id` col contain same id, then the reads of each sample will be merged.

+ For the same `assembly_group`, the reads will be combined together to do co-assembly,
  then the reads of each sample will map to the co-assembled results to calculate co-abundance for each `assembly_group`,
  then using MetaBAT2, MaxBin2, CONCOCT (You can specific which binner to use) to do binning.

+ For the same `binning_group`, the assembled results will be combined together to build minimap2 index,
  then the reads of each sample will map to the combined index to calculate co-abundance,
  then using vamb to do binning

+ `assembly_group` control the assembly and binning strategy; `binning_group` control the multisplit binning strategy developed by vamb.

+ If all `assembly_group` are unique, for MetaBAT2, MaxBin2, CONCOCT, etc., that is the traditional single-assembly and single-binning strategy.

+ If `binning_group` has multiple groups, then vamb will perform multisplit binning in each `binning_group` separately, which is useful for very large-scale data.
  For example, if there are 10,000 samples, and multisplit binning is performed once for every thousand samples, the binning group can be set to 10 different groups.

## Input format
- begin from trimming, rmhost or assembly:

  - `Paired-end reads`

  |  sample_id   |  assembly_group | binning_group |    fq1     |    fq2     |
  | :----------: | :-------------: | :-----------: | :--------: | :--------: |
  |  s1          | ag1             | bg1           | s1.1.fq.gz | s1.2.fq.gz |
  |  s2          | ag1             | bg1           | s2.1.fq.gz | s2.2.fq.gz |
  |  s3          | ag2             | bg2           | s3.1.fq.gz | s3.2.fq.gz |
  |  s4          | ag2             | bg2           | s4.1.fq.gz | s4.2.fq.gz |

  - `Paired-end reads(interleaved)`, update `config.yaml::reads_interleaved=true`  

  |  sample_id   |  assembly_group | binning_group |    fq1     |    fq2     |
  | :----------: | :-------------: | :-----------: | :--------: | :--------: |
  |  s1          | ag1             | bg1           | s1.fq.gz   |            |
  |  s2          | ag1             | bg1           | s2.fq.gz   |            |
  |  s3          | ag2             | bg2           | s3.fq.gz   |            |
  |  s4          | ag2             | bg2           | s4.fq.gz   |            |

  - (WIP) `Paired-end reads with long reads`, update `config.yaml::have_long=true`

  |  sample_id   |  assembly_group | binning_group |    fq1     |    fq2     |   fq_long      |
  | :----------: | :-------------: | :-----------: | :--------: | :--------: | :------------: |
  |  s1          | ag1             | bg1           | s1.1.fq.gz | s1.2.fq.gz |  s1.long.fq.gz |
  |  s2          | ag1             | bg1           | s2.1.fq.gz | s2.2.fq.gz |  s2.long.fq.gz |
  |  s3          | ag2             | bg2           | s3.1.fq.gz | s3.2.fq.gz |  s3.long.fq.gz |
  |  s4          | ag2             | bg2           | s4.1.fq.gz | s4.2.fq.gz |  s4.long.fq.gz |

  - (WIP)`Paired-end reads(interleaved) with long reads`, update `config.yaml::reads_interleaved=true` and `config.yaml::have_long=true`

  |  sample_id   |  assembly_group | binning_group |    fq1     |    fq2     |   fq_long      |
  | :----------: | :-------------: | :-----------: | :--------: | :--------: | :------------: |
  |  s1          | ag1             | bg1           | s1.fq.gz   |            |  s1.long.fq.gz |
  |  s2          | ag1             | bg1           | s2.fq.gz   |            |  s2.long.fq.gz |
  |  s3          | ag2             | bg2           | s3.fq.gz   |            |  s3.long.fq.gz |
  |  s4          | ag2             | bg2           | s4.fq.gz   |            |  s4.long.fq.gz |

  - `Single-end reads`, update `config.yaml::reads_layout=se`

  |  sample_id   |  assembly_group | binning_group |    fq1     |    fq2     |
  | :----------: | :-------------: | :-----------: | :--------: | :--------: |
  |  s1          | ag1             | bg1           | s1.fq.gz   |            |
  |  s2          | ag1             | bg1           | s2.fq.gz   |            |
  |  s3          | ag2             | bg2           | s3.fq.gz   |            |
  |  s4          | ag2             | bg2           | s4.fq.gz   |            |

  - `SRA (only support paired-end reads)` :
  SRA can be dumpped to Paired-end fastq reads

  |  sample_id   |  assembly_group | binning_group |    sra     |
  | :----------: | :-------------: | :-----------: | :--------: |
  |  s1          | ag1             | bg1           | s1.sra     |
  |  s2          | ag1             | bg1           | s2.sra     |
  |  s3          | ag2             | bg2           | s3.sra     |
  |  s4          | ag2             | bg2           | s4.sra     |

- begin from simulate, only support paired-end reads

  |  id   | genome | abundance | reads_num | model |
  | :---: | :----: | :-------: | :-------: | :---: |
  |  s1   | g1.fa  |    1.0    |    1M     | hiseq |
  |  s2   | g1.fa  |    0.5    |    2M     | hiseq |
  |  s2   | g2.fa  |    0.5    |    2M     | hiseq |
  |  s3   | g1.fa  |    0.2    |    3M     | hiseq |
  |  s3   | g2.fa  |    0.3    |    3M     | hiseq |
  |  s3   | g3.fa  |    0.5    |    3M     | hiseq |

It means:

The sample s1 contain 1M reads which come from g1, the relative abundance of
species g1 is 1.0.

The sample s2 contain 2M reads, 1M reads come from g1
and 1M reads come from g2. the relative abundance of
species g1 is 0.5, the relative abundance of
species g2 is 0.5.

The sample s3 contain 3M reads, 0.6M reads come from g1, 0.9M reads come from
g2 and 1.5M reads come from g3, the relative abundance of
species g1 is 0.2, the relative abundance of
species g2 is 0.3, the relative abundance of
species g3 is 0.5.

Then metapi will use [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) to generate metagenomics shotgun reads.

- (WIP) begin from binning
  |  sample_id   |  assembly_group | binning_group |    fq1     |    fq2     |   scaftigs          |
  | :----------: | :-------------: | :-----------: | :--------: | :--------: | :-----------------: |
  |  s1          | ag1             | bg1           | s1.1.fq.gz | s1.2.fq.gz |  ag1.scaftigs.fa.gz |
  |  s2          | ag1             | bg1           | s2.1.fq.gz | s2.2.fq.gz |  ag1.scaftigs.fa.gz |
  |  s3          | ag2             | bg2           | s3.1.fq.gz | s3.2.fq.gz |  ag2.scaftigs.fa.gz |
  |  s4          | ag2             | bg2           | s4.1.fq.gz | s4.2.fq.gz |  ag2.scaftigs.fa.gz |

- (WIP) begin from checkm, or dereplicate, or taxonomic assignment (fq1 and fq2 can be empty)
  |  sample_id   |   bins          |
  | :----------: | :-------------: |
  |  s1          |   s1.bin.1.fa   |
  |  s1          |   s1.bin.2.fa   |
  |  s1          |   s1.bin.3.fa   |
  |  s2          |   s2.bin.1.fa   |
  |  s2          |   s2.bin.2.fa   |
  |  s2          |   s3.bin.3.fa   |

## Default output structure (begin from trimming)
```
- config.yaml
- logs/
- profiles/
- envs/
- results/
    00.raw/
    01.trimming/
    02.rmhost/
    03.qcreport/
    04.assembly/
    05.alignment/
    06.binning/
    06.identify/
    07.predict/
    08.checkm/
    09.dereplicate/
    10.taxonomic/
    99.upload/
```

- We will try our best to keep the directory structure uniform. Sequence files are generally placed in the reads directory, and report files are generally placed in the report directory. 
- If you are not very clear about the output of the whole process, it is recommended to directly raise an issue or look at the code. 

## Getting help

If you want to report a bug or issue, or have problems with installing or
running the software, please create [a new
issue](https://github.com/ohmeta/metapi/issues).
This is the preferred way of getting support. Alternatively, you can [mail me](mailto:alienchuj@gmail.com).

## Contributing

Contributions welcome! Send me a pull request or get in [touch](mailto:alienchuj@gmail.com).

When contributing a PR, please use the [dev](https://github.com/ohmeta/metapi/tree/dev) branch.
For style, code will be checked using flake8,
[black](https://github.com/psf/black), and
[snakefmt](https://github.com/snakemake/snakefmt). These modules can be
installed via conda, `conda install black flake8 flake8-bugbear snakefmt` or via
pip `pip install black flake8 flake8-bugbear snakefmt`.

## Contributors

- Jie Zhu - [@alienzj](https://github.com/alienzj)
- Fangming Yang - [@yangfangming](https://github.com/yangfangming)
- Yanmei Ju - [@juyanmei](https://github.com/juyanmei)
- Weiting Liang - [@weiting-liang](https://github.com/weiting-liang)

## Citation

Over 50,000 Metagenomically Assembled Draft Genomes for the Human Oral Microbiome Reveal New Taxa  
Genomics, Proteomics & Bioinformatics, 2021, https://doi.org/10.1016/j.gpb.2021.05.001  

## License

This module is licensed under the terms of the [GPLv3 license](https://opensource.org/licenses/GPL-3.0).
