#!/usr/bin/env snakemake

import sys
import metapi
import pandas as pd
from snakemake.utils import min_version

min_version(7.0)

shell.executable("bash")

METAPI_DIR = metapi.__path__[0]
WRAPPER_DIR = os.path.join(METAPI_DIR, "wrappers")


RMHOST_DO = any([
    config["params"]["rmhost"]["bwa"]["do"],
    config["params"]["rmhost"]["bowtie2"]["do"]])


TRIMMING_DO = any([
    config["params"]["trimming"]["sickle"]["do"],
    config["params"]["trimming"]["fastp"]["do"],
    config["params"]["trimming"]["trimmomatic"]["do"]])


ASSEMBLERS = []
if config["params"]["assembly"]["megahit"]["do"]:
    ASSEMBLERS += ["megahit"]
if config["params"]["assembly"]["idba_ud"]["do"]:
    ASSEMBLERS += ["idba_ud"]
if config["params"]["assembly"]["metaspades"]["do"]:
    ASSEMBLERS += ["metaspades"]
if config["params"]["assembly"]["spades"]["do"]:
    ASSEMBLERS += ["spades"]


SAMPLES, DATA_TYPE = metapi.parse_samples(config["params"]["samples"])


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
        rules.raw_all.input,
        rules.trimming_all.input,
        rules.rmhost_all.input,
        rules.qcreport_all.input,
        rules.assembly_all.input,
        rules.predict_scaftigs_gene_all.input,
        rules.dereplicate_gene_all.input,
        rules.upload_all.input
