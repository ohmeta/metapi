<div align=center><img width="500" height="280" src="https://raw.githubusercontent.com/yangfangming/metapi/dev/docs/logo.svg"/></div>

# metapi

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![ohmeta-badge](https://img.shields.io/badge/install%20with-ohmeta-brightgreen.svg?style=flat)](http://anaconda.org/ohmeta)
[![PyPI version](https://badge.fury.io/py/metapi.svg)](https://badge.fury.io/py/metapi)
[![star this repo](http://githubbadges.com/star.svg?user=ohmeta&repo=metapi&style=flat)](https://github.com/ohmeta/metapi)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/metapi/badges/downloads.svg)](https://anaconda.org/bioconda/metapi)

A general metagenomics data mining system focus on robust microbiome research.

## Installation

metapi works with Python 3.6+.
You can install it via [bioconda](https://bioconda.github.io/):

```
# [WIP]
$ conda install -c bioconda metapi
# or [Complete]
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
```

### init

```
$ metapi init --help

  usage: metapi init [-h] [-d WORKDIR] [-s SAMPLES]
                    [-b {simulate,trimmingrmhost,assembly}]

  arguments:
      -h, --help            show this help message and exit
      -d, --workdir WORKDIR
                            project workdir, default: ./ (default: ./)
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

### Example

```
# init project
$ metapi init -d . -s samples.tsv -b trimming

# create conda environments (need connect to internet)
$ metapi mag_wf --conda_create_envs_only

# run pipeline with conda
# metapi mag_wf all --use_conda

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
$ metapi mag_wf profiling_metaphlan2_all --run --use_conda

# run MetaPhlAn3 profiling
$ metapi mag_wf profiling_metaphlan3_all --run

# run MAGs jgi profling (using jgi_summarize_bam_contig_depths)
$ metapi mag_wf profiling_jgi_all --run

# run HUMAnN2 profiling
$ metapi mag_wf profiling_humann2_all --run --use_conda

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

  | id  |    fq1     |    fq2     |
  | :-: | :--------: | :--------: |
  | s1  | aa.1.fq.gz | aa.2.fq.gz |
  | s2  | bb.1.fq.gz | bb.2.fq.gz |
  | s2  | cc.1.fq.gz | cc.2.fq.gz |
  | s3  | dd.1.fq.gz | dd.2.fq.gz |

  - `Paired-end reads(interleaved)`

  | id  |     fq1     | fq2 |
  | :-: | :---------: | :-: |
  | s1  | aa.12.fq.gz |     |
  | s2  | bb.12.fq.gz |     |
  | s2  | cc.12.fq.gz |     |
  | s3  | dd.12.fq.gz |     |

* `Paired-end reads with long reads`

| id  |    fq1     |    fq2     |    fq_long    |
| :-: | :--------: | :--------: | :-----------: |
| s1  | aa.1.fq.gz | aa.2.fq.gz | aa.long.fq.gz |
| s2  | bb.1.fq.gz | bb.2.fq.gz | bb.long.fq.gz |
| s2  | cc.1.fq.gz | cc.2.fq.gz | cc.long.fq.gz |
| s3  | dd.1.fq.gz | dd.2.fq.gz | dd.long.fq.gz |

- `Paired-end reads(interleaved) with long reads`

| id  |     fq1     | fq2 |    fq_long    |
| :-: | :---------: | :-: | :-----------: |
| s1  | aa.12.fq.gz |     | aa.long.fq.gz |
| s2  | bb.12.fq.gz |     | bb.long.fq.gz |
| s2  | cc.12.fq.gz |     | cc.long.fq.gz |
| s3  | dd.12.fq.gz |     | dd.long.fq.gz |

- `Single-end reads`

| id  |    fq1     | fq2 |
| :-: | :--------: | :-: |
| s1  | aa.1.fq.gz |     |
| s2  | bb.1.fq.gz |     |
| s2  | cc.1.fq.gz |     |
| s3  | dd.1.fq.gz |     |

- `SRA (only support paired-end reads)` :

SRA can be dumpped to Paired-end fastq reads

| id  |  sra   |
| :-: | :----: |
| s1  | aa.sra |
| s2  | bb.sra |
| s2  | cc.sra |
| s3  | dd.sra |

- begin from simulate, only support paired-end reads

  | id  | genome | abundance | reads_num | model |
  | :-: | :----: | :-------: | :-------: | :---: |
  | s1  | g1.fa  |    1.0    |    1M     | hiseq |
  | s2  | g1.fa  |    0.5    |    2M     | hiseq |
  | s2  | g2.fa  |    0.5    |    2M     | hiseq |
  | s3  | g1.fa  |    0.2    |    3M     | hiseq |
  | s3  | g2.fa  |    0.3    |    3M     | hiseq |
  | s3  | g3.fa  |    0.5    |    3M     | hiseq |

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
