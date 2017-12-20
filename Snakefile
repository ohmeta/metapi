#!/usr/bin/env snakemake
configfile: "config/config.yaml"
import os
import json
from scripts.find_path import find_path_tag

(R1_raw, R2_raw) = find_path_tag(config["raw_reads_dir"], "raw")
assert sorted(R2_raw.keys()) == sorted(R2_raw.keys())
config["R1_raw"] = R1_raw
config["R2_raw"] = R2_raw
print(json.dumps(config, sort_keys=True, indent=4))

include: "rules/quality_control/filter_zebra.rules"
#include: "rules/quality_control/rmhost.rules"
#include: "rules/quality_control/qc_report.rules"
