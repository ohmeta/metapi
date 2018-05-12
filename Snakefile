#!/usr/bin/env snakemake
import pandas as pd
shell.executable("bash")

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

rule all:
    input:
        expand("{trim}/{unit.sample}_{unit.unit}.trimmed.{read}.fq.gz",
               trim=config["results"]["trim"],
               unit=units.reset_index().itertuples(),
               read=[1, 2, 'single'])

include: "rules/trim.smk"
