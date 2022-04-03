#!/usr/bin/env snakemake

import sys
import metapi
import pandas as pd
from pprint import pprint

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


ASSEMBLERS_CO = []
if config["params"]["coassembly"]["megahit"]["do"]:
    ASSEMBLERS_CO += ["megahit"]


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


include: "../rules/simulate.smk"
include: "../rules/raw.smk"
include: "../rules/trimming.smk"
include: "../rules/rmhost.smk"
include: "../rules/qcreport.smk"
include: "../rules/assembly.smk"
#include: "../rules/coassembly.smk"
#include: "../rules/predict_scaftigs.smk"
#include: "../rules/copredict_scaftigs.smk"
#include: "../rules/alignment.smk"
#include: "../rules/coalignment.smk"
#include: "../rules/binning.smk"
#include: "../rules/cobinning.smk"
#include: "../rules/multisplit_binning.smk"
#include: "../rules/predict_bins.smk"
#include: "../rules/copredict_bins.smk"
#include: "../rules/checkm.smk"
#include: "../rules/cocheckm.smk"
#include: "../rules/dereplicate_mags.smk"
#include: "../rules/classify.smk"
#include: "../rules/coclassify.smk"
#include: "../rules/profiling.smk"
#include: "../rules/upload.smk"


rule all:
    input:
        rules.simulate_all.input,
        rules.raw_all.input,
        rules.trimming_all.input,
        rules.rmhost_all.input,
        rules.qcreport_all.input,
        rules.assembly_all.input#,
#        rules.predict_scaftigs_gene_all.input,
#        rules.alignment_all.input,
#        rules.binning_all.input,
#        rules.predict_bins_gene_all.input,
#        rules.checkm_all.input,
#        rules.dereplicate_mags_all.input,
#        rules.classify_all.input,
#        rules.profiling_all.input,
#        rules.upload_all.input


localrules: \
    simulate_all, \
    prepare_short_reads_all, \
    prepare_long_reads_all, \
    prepare_reads_all, \
    raw_fastqc_all, \
    raw_report_all, \
    trimming_oas1_all, \
    trimming_sickle_all, \
    trimming_fastp_all, \
    trimming_report_all, \
    trimming_all, \
    rmhost_soap_all, \
    rmhost_bwa_all, \
    rmhost_bowtie2_all, \
    rmhost_minimap2_all, \
    rmhost_kraken2_all, \
    rmhost_report_all, \
    rmhost_all, \
    qcreport_all, \
    assembly_megahit_all, \
    assembly_idba_ud_all, \
    assembly_metaspades_all, \
    assembly_spades_all, \
    assembly_plass_all, \
    assembly_opera_ms_all, \
    assembly_metaquast_all, \
    assembly_report_all, \
    single_assembly_all, \
    coassembly_megahit_all, \
    coassembly_all, \
    assembly_all, \
    alignment_base_depth_all, \
    coalignment_base_depth_all, \
    alignment_report_all, \
    single_alignment_all, \
    coalignment_report_all, \
    coalignment_all, \
    alignment_all, \
    binning_metabat2_coverage_all, \
    binning_metabat2_all, \
    binning_maxbin2_all, \
    binning_concoct_all, \
    binning_graphbin2_all, \
    binning_dastools_all, \
    binning_vamb_prepare_all, \
    binning_vamb_all, \
    binning_report_all, \
    single_binning_all, \
    cobinning_metabat2_coverage_all, \
    cobinning_metabat2_all, \
    cobinning_maxbin2_all, \
    cobinning_concoct_all, \
    cobinning_graphbin2_all, \
    cobinning_dastools_all, \
    cobinning_report_all, \
    cobinning_all, \
    multisplit_binning_all, \
    binning_all, \
    single_predict_scaftigs_gene_all, \
    single_predict_bins_gene_all, \
    copredict_scaftigs_gene_all, \
    copredict_bins_gene_all, \
    predict_scaftigs_gene_all, \
    predict_bins_gene_all, \
    copredict_all, \
    predict_all, \
    single_checkm_all, \
    cocheckm_all, \
    checkm_all, \
    dereplicate_mags_drep_all, \
    dereplicate_mags_all, \
    dereplicate_all, \
    classify_short_reads_kraken2_krona_report,
    classify_short_reads_kraken2_combine_kreport,
    classify_short_reads_kraken2_combine_kreport_mpa,
    classify_short_reads_kraken2_all, \
    single_classify_hmq_bins_gtdbtk_all, \
    coclassify_hmq_bins_gtdbtk_all, \
    classify_hmq_bins_gtdbtk_all, \
    single_classify_all, \
    coclassify_all, \
    classify_all, \
    profiling_bgi_soap_all, \
    profiling_bowtie2_all, \
    profiling_metaphlan2_all, \
    profiling_metaphlan3_all, \
    profiling_jgi_all, \
    profiling_bracken_all, \
    profiling_humann2_config, \
    profiling_humann2_all, \
    profiling_humann3_config, \
    profiling_humann3_all, \
    profiling_all, \
    upload_sequencing_all, \
    upload_assembly_all, \
    upload_all
