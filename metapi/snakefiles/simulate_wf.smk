
#!/usr/bin/env snakemake

import sys
from pprint import pprint
import pandas as pd

from snakemake.utils import min_version
min_version("7.0")
shell.executable("bash")

import metapi

METAPI_DIR = metapi.__path__[0]
WRAPPER_DIR = os.path.join(METAPI_DIR, "wrappers")
DATA_DIR = os.path.join(METAPI_DIR, "data")


SAMPLES, DATA_TYPE = metapi.parse_samples(config["params"]["samples"])


include: "../rules/simulate.smk"


rule all:
    input:
        rules.simulate_all.input