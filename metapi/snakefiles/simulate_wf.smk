
#!/usr/bin/env snakemake

import sys
from pprint import pprint
import pandas as pd

#import metapi

from snakemake.utils import min_version
min_version("7.0")
shell.executable("bash")

sys.path.insert(0, "/home/jiezhu/toolkit/metapi_dev")
import metapi

METAPI_DIR = metapi.__path__[0]
WRAPPER_DIR = os.path.join(METAPI_DIR, "wrappers")
DATA_DIR = os.path.join(METAPI_DIR, "data")

pprint(METAPI_DIR)


SAMPLES = metapi.parse_genomes(config["params"]["samples"],
                               config["output"]["simulate"])


include: "../rules/simulate.smk"


rule all:
    input:
        rules.simulate_all.input