[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![ohmeta-badge](https://img.shields.io/badge/install%20with-ohmeta-brightgreen.svg?style=flat)](http://anaconda.org/ohmeta)
[![PyPI version](https://badge.fury.io/py/metapi.svg)](https://badge.fury.io/py/metapi)
[![star this repo](http://githubbadges.com/star.svg?user=ohmeta&repo=metapi&style=flat)](https://github.com/ohmeta/metapi)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metapi/badges/downloads.svg)](https://anaconda.org/bioconda/metapi)

<div align=center><img width="500" height="280" src="docs/logo.svg"/></div>

A general metagenomics data mining system focus on robust microbiome research.

## Overview
### MAG workflow
<div align=center><img width="600" height="800" src="docs/mag_workflow.svg"></div>

## Installation

metapi works with Python 3.6+.
You can install it via [bioconda](https://bioconda.github.io/):

```
$ conda install -c bioconda metapi
# or
$ conda install -c ohmeta metapi
```

Or via pip:

```
$ pip install metapi
```

## Run

### help

```
$ metapi --help

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

### init

```
$ metapi init --help

  usage: metapi init [-h] [-d WORKDIR] [-s SAMPLES]
                    [-b {simulate,trimmingrmhost,assembly}]

  arguments:
      -h, --help            show this help message and exit
      -d, --workdir WORKDIR
                            project workdir (default: ./)
      -s, --samples SAMPLES
                            desired input:
                            samples list, tsv format required.

                            if begin from trimming, rmhost, or assembly:
                                if it is fastq:
                                    the header is [id, fq1, fq2]
                                if it is sra:
                                    the header is [id, sra]
                            if begin from simulate:
                                the header is [id, genome, abundance, reads_num, model]

    -b, --begin {simulate,trimming,rmhost,assembly}
                            pipeline starting point (default: trimming)
```

### mag_wf

```
$ metapi mag_wf --help

  usage: metapi mag_wf [-h] [-d WORKDIR] [--config CONFIG] [--cluster CLUSTER]
                            [--cores CORES] [--jobs JOBS] [--list] [--run] [--debug]
                            [--dry-run] [--qsub] [--wait WAIT] [--use-conda]
                            [--conda-prefix CONDA_PREFIX] [--conda-create-envs-only]
                            [TASK]

  positional arguments:
  TASK              pipeline end point. Allowed values are 
       simulate_all, prepare_short_reads_all, prepare_long_reads_all,
       prepare_reads_all,
       raw_fastqc_all, raw_report_all, raw_all,
       trimming_oas1_all, trimming_sickle_all,
       trimming_fastp_all, trimming_report_all,
       trimming_all,
       rmhost_soap_all, rmhost_bwa_all, rmhost_bowtie2_all,
       rmhost_minimap2_all, rmhost_kraken2_all,
       rmhost_report_all, rmhost_all, qcreport_all,
       assebmly_megahit_all, assembly_idba_ud_all,
       assembly_metaspades_all, assembly_spades_all,
       assembly_plass_all, assembly_opera_ms_all,
       assembly_metaquast_all, assembly_report_all,
       single_assembly_all, coassembly_megahit_all,
       coassembly_all, assembly_all,
       alignment_base_depth_all, single_alignment_all,
       coalignment_base_depth_all,
       coalignment_all, alignment_all,
       binning_metabat2_coverage_all, binning_metabat2_all,
       binning_maxbin2_all, binning_concoct_all, binning_graphbin2_all,
       binning_dastools_all, binning_vamb_all,
       binning_report_all, single_binning_all,
       cobinning_metabat2_coverage_all, cobinning_metabat2_all,
       cobinning_maxbin2_all, cobinning_concoct_all,
       cobinning_graphbin2_all, cobinning_dastools_all, 
       cobinning_report_all, cobinning_all, binning_all,
       predict_scaftigs_gene_prodigal_all, predict_scaftigs_gene_prokka_all,
       predict_bins_gene_prodigal_all, predict_bins_gene_prokka_all,
       single_predict_scaftigs_gene_all, single_predict_bins_gene_all,
       copredict_scaftigs_gene_prodigal_all,
       copredict_scaftigs_gene_prokka_all,
       copredict_bins_gene_prodigal_all, copredict_bins_gene_prokka_all,
       copredict_scafitgs_gene_all, copredict_bins_gene_all,
       predict_scaftigs_gene_all, predict_bins_gene_all,
       copredict_all, predict_all,
       single_checkm_all, cocheckm_all, checkm_all,
       dereplicate_mags_drep_all, dereplicate_mags_all,
       classify_short_reads_kraken2_all,
       single_classify_hmq_bins_gtdbtk_all,
       coclassify_hmq_bins_gtdbtk_all, classify_hmq_bins_gtdbtk_all,
       single_classify_all, coclassify_all, classify_all,
       profiling_bgi_soap_all, profiling_bowtie2_all,
       profiling_metaphlan2_all, profiling_metaphlan3_all,
       profiling_jgi_all, profiling_bracken_all,
       profiling_humann2_all, profiling_humann3_all,
       profiling_all, upload_sequencing_all, upload_assembly_all, upload_all,
       all (default: all)

  optional arguments:
  -h, --help            show this help message and exit
  -d, --workdir WORKDIR
                        project workdir (default: ./)
  --config CONFIG       config.yaml (default: ./config.yaml)
  --cluster CLUSTER     cluster.yaml (default: ./cluster.yaml)
  --cores CORES         CPU cores (default: 8)
  --jobs JOBS           qsub job numbers (default: 80)
  --list                list pipeline rules
  --run                 run pipeline
  --debug               debug pipeline
  --dry-run             dry run pipeline
  --qsub                qsub pipeline
  --wait WAIT           wait given seconds (default: 60)
  --use-conda           use conda environment
  --conda-prefix CONDA_PREFIX
                        conda environment prefix
                        (default: /ldfssz1/ST_META/share/User/zhujie/.conda/envs)
  --conda-create-envs-only
    conda create environments only
```

### Example

```
# init project
$ metapi init -d . -s samples.tsv -b trimming

# create conda environments
$ metapi mag_wf --conda-create-envs-only

# run pipeline with conda
# metapi mag_wf all --use-conda
# metapi mag_wf all --use-conda --conda-prefix /path/to/your/default/envs/dir

# run raw_fastqc
$ metapi mag_wf raw_fastqc_all --run

# run trimming
$ metapi mag_wf trimming_all --run

# run rmhost
$ metapi mag_wf rmhost_all --run

# run qc report
$ metapi mag_wf qcreport_all --run

# run assembly
$ metapi mag_wf assembly_all --run

# run binning
$ metapi mag_wf binning_all --run

# run gene predict
$ metapi mag_wf predict_all --run

# run MAGs checkm
$ metapi mag_wf checkm_all --run

# run MAGs classify
$ metapi mag_wf classify_all --run

# run MetaPhlAn2 profiling
$ metapi mag_wf profiling_metaphlan2_all --run \
  --use-conda --conda-prefix /ldfssz1/ST_META/share/User/zhujie/.conda/envs

# run MetaPhlAn3 profiling
$ metapi mag_wf profiling_metaphlan3_all --run

# run MAGs jgi profling (using jgi_summarize_bam_contig_depths)
$ metapi mag_wf profiling_jgi_all --run

# run HUMAnN2 profiling
$ metapi mag_wf profiling_humann2_all --run \
  --use-conda --conda-prefix /ldfssz1/ST_META/share/User/zhujie/.conda/envs

# run mag_wf all
$ metapi mag_wf --run

# run gene_wf all
$ metapi gene_wf --run
```

## input requirements

The input samples file: `samples.tsv` format:

Note: If `id` col contain same id, then the reads of each sample will be merged.
Note: The fastq need gzip compress.

- begin from trimming, rmhost or assembly:

  - `Paired-end reads`

  |  id   |    fq1     |    fq2     |
  | :---: | :--------: | :--------: |
  |  s1   | aa.1.fq.gz | aa.2.fq.gz |
  |  s2   | bb.1.fq.gz | bb.2.fq.gz |
  |  s2   | cc.1.fq.gz | cc.2.fq.gz |
  |  s3   | dd.1.fq.gz | dd.2.fq.gz |

  - `Paired-end reads(interleaved)`

  |  id   |     fq1     |  fq2  |
  | :---: | :---------: | :---: |
  |  s1   | aa.12.fq.gz |       |
  |  s2   | bb.12.fq.gz |       |
  |  s2   | cc.12.fq.gz |       |
  |  s3   | dd.12.fq.gz |       |

* `Paired-end reads with long reads`

|  id   |    fq1     |    fq2     |    fq_long    |
| :---: | :--------: | :--------: | :-----------: |
|  s1   | aa.1.fq.gz | aa.2.fq.gz | aa.long.fq.gz |
|  s2   | bb.1.fq.gz | bb.2.fq.gz | bb.long.fq.gz |
|  s2   | cc.1.fq.gz | cc.2.fq.gz | cc.long.fq.gz |
|  s3   | dd.1.fq.gz | dd.2.fq.gz | dd.long.fq.gz |

- `Paired-end reads(interleaved) with long reads`

|  id   |     fq1     |  fq2  |    fq_long    |
| :---: | :---------: | :---: | :-----------: |
|  s1   | aa.12.fq.gz |       | aa.long.fq.gz |
|  s2   | bb.12.fq.gz |       | bb.long.fq.gz |
|  s2   | cc.12.fq.gz |       | cc.long.fq.gz |
|  s3   | dd.12.fq.gz |       | dd.long.fq.gz |

- `Single-end reads`

|  id   |    fq1     |  fq2  |
| :---: | :--------: | :---: |
|  s1   | aa.1.fq.gz |       |
|  s2   | bb.1.fq.gz |       |
|  s2   | cc.1.fq.gz |       |
|  s3   | dd.1.fq.gz |       |

- `SRA (only support paired-end reads)` :

SRA can be dumpped to Paired-end fastq reads

|  id   |  sra   |
| :---: | :----: |
|  s1   | aa.sra |
|  s2   | bb.sra |
|  s2   | cc.sra |
|  s3   | dd.sra |

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

The sample s1 contain 1M reads which come from g1, the relatative abundance of
species g1 is 1.0.

The sample s2 contain 2M reads, 1M reads come from g1
and 1M reads come from g2. the relatative abundance of
species g1 is 0.5, the relatative abundance of
species g2 is 0.5.

The sample s3 contain 3M reads, 0.6M reads come from g1, 0.9M reads come from
g2 and 1.5M reads come from g3, the relatative abundance of
species g1 is 0.2, the relatative abundance of
species g2 is 0.3, the relatative abundance of
species g3 is 0.5.

Then metapi will use [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) to generate metagenomics shutgun reads.

## FAQ

- You know what you want to do, so you know how to configure config.yaml
- You know snakemake, so you know how to hack metapi

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

## License

This module is licensed under the terms of the [GPLv3 license](https://opensource.org/licenses/GPL-3.0).
