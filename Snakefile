#!/usr/bin/env snakemake
configfile: "config.json"
import os

file_ending = config["reads_file_ending"]
R1_reads = {}
R2_reads = {}

for reads_f in os.listdir(config["raw_reads_dir"]):
    if reads_f.endswith("1." + file_ending):
        sample_name = reads_f.rstrip(".|_|-" + "1." + file_ending)
        R1_reads[sample_name] = os.path.join(config["raw_reads_dir"], reads_f)
    elif reads_f.endswith("2." + file_ending):
        sample_name = reads_f.rstrip(".|_|-" + "2." + file_ending)
        R2_reads[sample_name] = os.path.join(config["raw_reads_dir"], reads_f)

assert sorted(R1_reads.keys()) == sorted(R2_reads.keys())

config["R1_reads"] = sorted(R1_reads.values())
config["R2_reads"] = sorted(R2_reads.values())
config["R1_reads_per_sample"] = R1_reads
config["R2_reads_per_sample"] = R2_reads

include: "rules/quality_control/filter_zebra.rules"
