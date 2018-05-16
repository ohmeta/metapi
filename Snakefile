#!/usr/bin/env snakemake
import os
import sys
import shutil
import pandas as pd

shell.executable("bash")

program_list = ["sickle", "bwa", "samtools", "megahit", "metabat2", "pigz"]

def checking_dependencies(program_list):
    install = []
    exit = False
    for program in program_list:
        where = shutil.which(program)
        if where is None:
            exit = True
            install.append(program)
            print(program + ":\tno")
        else:
            print(program + ":\tyes")

    if exit:
        if "metabat2" not in install:
            print("\npelase use conda to install these program:")
            print("conda install %s" % (" ".join(install)))
        else:
            install_info = " ".join(install).replace("metabat2", "")
            print("\npelase use conda to install these program:")
            if install_info != "":
                print("conda install %s" % install_info)
            print("conda install -c ursky metabat2")
        sys.exit()

checking_dependencies(program_list)

configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col=["sample"])
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])


def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def debug(samples, units):
    print(samples)
    print(samples.index)
    print("\n")
    print(units)
    print(units.index)
    print("\n")
    for unit in units.reset_index().itertuples():
        print(unit)
        print(unit.sample)
        print(unit.unit)

#debug(samples, units)

# test trim
'''
rule all:
    input:
        expand("{trim}/{unit.sample}_{unit.unit}.trimmed.{read}.fq.gz",
               trim=config["results"]["trim"],
               unit=units.reset_index().itertuples(),
               read=[1, 2, 'single'])
'''

# test rmhost
'''
rule all:
    input:
        expand(["{rmhost}/{unit.sample}_{unit.unit}.rmhost.{read}.fq.gz",
                "{rmhost}/{unit.sample}_{unit.unit}.flagstat.txt"],
               rmhost=config["results"]["rmhost"],
               unit=units.reset_index().itertuples(),
               read=["1", "2"])
'''

# test individual assembly
'''
rule all:
    input:
        expand("{assembly}/{unit.sample}_{unit.unit}.megahit_out/{unit.sample}_{unit.unit}.contigs.fa",
               assembly=config["results"]["assembly"],
               unit=units.reset_index().itertuples())
'''

# test algnment
'''
rule all:
    input:
        expand(["{alignment}/{unit.sample}_{unit.unit}.sorted.bam",
                "{alignment}/{unit.sample}_{unit.unit}.flagstat.txt"],
               alignment=config["results"]["alignment"],
               unit=units.reset_index().itertuples())
'''

# test binning
rule all:
    input:
        expand("{binning}/bins/{unit.sample}_{unit.unit}.metabat2_out/done",
               binning=config["results"]["binning"],
               unit=units.reset_index().itertuples())


include: "rules/trim.smk"
include: "rules/rmhost.smk"
#include: "rules/qcreport.smk"
include: "rules/assembly.smk"
include: "rules/alignment.smk"
include: "rules/binning.smk"
