#!/usr/bin/env snakemake

import sys
import metapi
import pandas as pd
from pprint import pprint
from snakemake.utils import min_version

min_version(7.0)

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


TRIMMING_DO = True \
    if config["params"]["trimming"]["oas1"]["do"] or \
       config["params"]["trimming"]["sickle"]["do"] or \
       config["params"]["trimming"]["fastp"]["do"] \
       else False


RMHOST_DO = True \
    if config["params"]["rmhost"]["soap"]["do"] or \ 
       config["params"]["rmhost"]["bwa"]["do"] or \
       config["params"]["rmhost"]["bowtie2"]["do"] or \
       config["params"]["rmhost"]["minimap2"]["do"] or \
       config["params"]["rmhost"]["kraken2"]["do"] or \
       config["params"]["rmhost"]["kneaddata"]["do"] \
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

if config["params"]["assembly"]["opera_ms"]["do"]:
    ASSEMBLERS += ["opera_ms"]
    if (config["params"]["assembly"]["opera_ms"]["short_read_assembler"] == "megahit") \
       and (not "megahit" in ASSEMBLERS):
        config["params"]["assembly"]["megahit"]["do"] = True
        ASSEMBLERS += ["megahit"]
    elif (config["params"]["assembly"]["opera_ms"]["short_read_assembler"] == "metaspades") \
         and (not "metaspades" in ASSEMBLERS):
        config["params"]["assembly"]["metaspades"]["do"] = True
        ASSEMBLERS += ["metaspades"]


BINNERS_TOTAL = []
BINNERS_GRAPHBIN = []
BINNERS_DASTOOLS = []

if config["params"]["binning"]["metabat2"]["do"]:
    BINNERS_TOTAL += ["metabat2"]
if config["params"]["binning"]["maxbin2"]["do"]:
    BINNERS_TOTAL += ["maxbin2"]
if config["params"]["binning"]["concoct"]["do"]:
    BINNERS_TOTAL += ["concoct"]
if config["params"]["binning"]["vamb"]["do"]:
    BINNERS_TOTAL += ["vamb"]

if config["params"]["binning"]["graphbin2"]["do"]:
    BINNERS_GRAPHBIN = BINNERS_TOTAL.copy()
    for i in BINNERS_GRAPHBIN:
        BINNERS_TOTAL.append(i + "_graphbin2")
        BINNERS_DASTOOLS.append(i + "_graphbin2")
else:
    BINNERS_DASTOOLS = BINNERS_TOTAL.copy()

if config["params"]["binning"]["dastools"]["do"]:
    BINNERS_TOTAL.append("dastools")


BINNERS_CHECKM = config["params"]["checkm"]["check_binners"]


if config["params"]["simulate"]["do"]:
    SAMPLES = metapi.parse_genomes(config["params"]["samples"],
                                   config["output"]["simulate"])
else:
    SAMPLES = metapi.parse_samples(config["params"]["samples"],
                                   config["params"]["interleaved"],
                                   config["params"]["reads_layout"],
                                   config["params"]["begin"])

    SAMPLES_ID_LIST = SAMPLES.index.get_level_values("sample_id").unique()
    SAMPLES_ASSEMBLY_GROUP_LIST = SAMPLES.index.get_level_values("assembly_group").unique()
    SAMPLES_BINNING_GROUP_LIST = SAMPLES.index.get_level_values("binning_group").unique()


READS_FORMAT = "sra" \
    if "sra" in SAMPLES.columns \
       else "fastq"


## TODO
"""
if config["params"]["begin"] == "binning":
    for sample_id in SAMPLES.index.unique():
        scaftigs = os.path.abspath(SAMPLES.loc[sample_id, "scaftigs"])
        if len(ASSEMBLERS) == 0:
            print("When begin with binning, please specific assembler")
        else:
            for assembler in ASSEMBLERS:
                scaftigs_ = os.path.join(
                    os.path.abspath(config["output"]["assembly"]),
                    f"scaftigs/{sample_id}.{assembler}.out/{sample_id}.{assembler}.scaftigs.fa.gz")
                if not os.path.exists(scaftigs_):
                    scaftigs_dir = os.path.dirname(scaftigs_)
                    shell(f'''mkdir -p {scaftigs_dir}''')
                    shell(f'''ln -s {scaftigs} {scaftigs_}''')
"""


include: "../rules/simulate.smk"
include: "../rules/raw.smk"
include: "../rules/trimming.smk"
include: "../rules/rmhost.smk"
include: "../rules/qcreport.smk"
include: "../rules/assembly.smk"
include: "../rules/predict_scaftigs.smk"
include: "../rules/alignment.smk"
include: "../rules/binning.smk"
include: "../rules/binning_multisplit.smk"
include: "../rules/binning_refine.smk"
include: "../rules/binning_report.smk"
include: "../rules/identify_single.smk"
include: "../rules/identify_multi.smk"
include: "../rules/predict_bins.smk"
include: "../rules/checkm.smk"
include: "../rules/dereplicate_mags.smk"
include: "../rules/taxonomic.smk"
include: "../rules/upload.smk"


rule all:
    input:
        rules.simulate_all.input,
        rules.raw_all.input,
        rules.trimming_all.input,
        rules.rmhost_all.input,
        rules.qcreport_all.input,
        rules.assembly_all.input,
        rules.predict_scaftigs_gene_all.input,
        rules.alignment_all.input,
        rules.binning_all.input,
        rules.identify_all.input,
        rules.predict_bins_gene_all.input,
        rules.checkm_all.input,
        rules.dereplicate_mags_all.input,
        rules.taxonomic_all.input,
        rules.upload_all.input