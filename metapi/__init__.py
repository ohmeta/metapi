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
from metapi.sampler import get_samples_id_by_assembly_group
from metapi.sampler import get_samples_id_by_binning_group
from metapi.sampler import get_assembly_group_by_binning_group
from metapi.sampler import get_binning_group_by_assembly_group
from metapi.sampler import get_multibinning_group_by_assembly_group

from metapi.qcer import change
from metapi.qcer import compute_host_rate
from metapi.qcer import qc_bar_plot
from metapi.qcer import parse_fastp_json

from metapi.assembler import assembler_init
from metapi.assembler import parse_assembly

from metapi.aligner import flagstats_summary

from metapi.predictor import parse_gff
from metapi.predictor import extract_faa

from metapi.binner import get_binning_info
from metapi.binner import generate_bins
from metapi.binner import extract_bins_report
from metapi.binner import combine_jgi

from metapi.checkmer import checkm_prepare
from metapi.checkmer import checkm_reporter

from metapi.classifier import demultiplex
from metapi.classifier import gtdbtk_prepare

from metapi.uploader import gen_samples_info
from metapi.uploader import gen_info

from metapi.__about__ import __version__, __author__

name = "metapi"
