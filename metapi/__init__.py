#!/usr/bin/env python

from metapi.configer import metaconfig
from metapi.configer import parse_yaml
from metapi.configer import update_config
from metapi.configer import custom_help_formatter

from metapi.tooler import parse
from metapi.tooler import merge

from metapi.simulator import parse_genomes
from metapi.simulator import get_simulate_info
from metapi.simulator import simulate_short_reads

from metapi.sampler import parse_samples
from metapi.sampler import get_reads
from metapi.sampler import get_sample_id
from metapi.sampler import get_sample_id_

from metapi.qcer import change
from metapi.qcer import compute_host_rate
from metapi.qcer import qc_bar_plot

from metapi.assembler import assembler_init
from metapi.assembler import parse_assembly

from metapi.aligner import flagstats_summary

from metapi.binner import get_binning_info
from metapi.binner import generate_bins

from metapi.checkmer import checkm_report

from metapi.classifier import demultiplex

from metapi.profiler import profiler_init
from metapi.profiler import get_all_abun_df
from metapi.profiler import get_profile
from metapi.profiler import metaphlan_init
from metapi.profiler import merge_metaphlan_tables

from metapi.uploader import gen_samples_info
from metapi.uploader import gen_info

from metapi.__about__ import __version__, __author__

name = "metapi"
