#!/usr/bin/env snakemake
import pandas as pd
shell.executable("bash")

configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i inunits.index.levels])

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit)], "2.fq.gz")

include: "rules/trim.smk"
include: "rules/assembly.smk"
include: "rules/align.smk"
include: "rules/binning.smk"
