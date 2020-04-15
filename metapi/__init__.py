#!/usr/bin/env python

from metapi.tooler import parse
from metapi.tooler import merge

from metapi.sampler import parse_samples
from metapi.sampler import parse_cobin_samples_id
from metapi.sampler import get_reads
from metapi.sampler import get_sample_id

from metapi.qcer import change
from metapi.qcer import compute_host_rate

from metapi.assembler import assembler_init
from metapi.assembler import parse_assembly

from metapi.aligner import flagstats_summary

from metapi.checkmer import report

from metapi.classifier import demultiplex

from metapi.profiler import profiler_init
from metapi.profiler import get_all_abun_df
from metapi.profiler import get_profile

from metapi.uploader import gen_samples_info
from metapi.uploader import gen_info

name = "metapi"
__version__ = "0.7.0"
