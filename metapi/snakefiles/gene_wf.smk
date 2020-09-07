#!/usr/bin/env snakemake

import sys
import metapi
import pandas as pd

shell.executable("bash")

METAPI_DIR = metapi.__path__[0]
WRAPPER_DIR = os.path.join(METAPI_DIR, "wrappers")


IS_PE = True \
    if config["params"]["reads_layout"] == "pe" \
       else False


IS_INTERLEAVED = True \
    if config["params"]["interleaved"] \
       else False


HAVE_LONG = True \
    if IS_PE and config["params"]["have_long"] \
       else False


RMHOST_DO = True \
    if config["params"]["rmhost"]["bwa"]["do"] or \
       config["params"]["rmhost"]["bowtie2"]["do"] \
       else False


TRIMMING_DO = True \
    if config["params"]["trimming"]["oas1"]["do"] or \
       config["params"]["trimming"]["sickle"]["do"] or \
       config["params"]["trimming"]["fastp"]["do"] \
       else False


ASSEMBLERS = []
if config["params"]["assembly"]["megahit"]["do"]:
    ASSEMBLERS += ["megahit"]
if config["params"]["assembly"]["idba_ud"]["do"]:
    ASSEMBLERS += ["idba_ud"]
if config["params"]["assembly"]["metaspades"]["do"]:
    ASSEMBLERS += ["metaspades"]
if config["params"]["assembly"]["spades"]["do"]:
    ASSEMBLERS += ["spades"]

if config["params"]["simulate"]["do"]:
    SAMPLES = metapi.parse_genomes(config)
else:
    SAMPLES = metapi.parse_samples(config)

READS_FORMAT = "sra" \
    if "sra" in SAMPLES.columns \
       else "fastq"


include: "../rules/simulate.smk"
include: "../rules/raw.smk"
include: "../rules/trimming.smk"
include: "../rules/rmhost.smk"
include: "../rules/qcreport.smk"
include: "../rules/assembly.smk"
include: "../rules/predict_scaftigs.smk"
include: "../rules/dereplicate_cds.smk"
include: "../rules/upload.smk"


rule all:
    input:
        rules.simulate_all.input,
        rules.raw_all.input,
        rules.trimming_all.input,
        rules.rmhost_all.input,
        rules.qcreport_all.input,
        rules.single_assembly_all.input,
        rules.single_predict_scaftigs_gene_all.input,
        rules.dereplicate_gene_all.input,
        rules.upload_all.input
