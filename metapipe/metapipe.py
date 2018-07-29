#!/usr/bin/env python

import argparse
import logging
import os
import shutil
import subprocess
import sys

import yaml

import metaconfig
import metasample
from snakemake.utils import update_config

__version__ = "0.1.0"

simulation_steps = [
    "simulation", "fastqc", "trim", "rmhost", "qc_report", "assembly",
    "alignment", "binning", "checkm", "dereplication", "classification",
    "annotation"
]

workflow_steps = [
    "fastqc", "trim", "rmhost", "qc_report", "assembly",
    "alignment", "binning", "checkm", "dereplication", "classification",
    "annotation"
]


def initialization(args):
    if args.workdir:
        proj = metaconfig.config(args.workdir)
        proj.create_dirs()
        default_config = proj.get_config()
        #print(type(default_config))

        #config = {}
        #config["params"] = {}
        #config["params"]["cluster"] = {}
        #config["results"] = {}
        #config["results"]["raw"] = {}

        if args.queue:
            default_config["params"]["cluster"]["queue"] = args.queue
        if args.project:
            default_config["params"]["cluster"]["project"] = args.project
        if args.samples:
            default_config["results"]["raw"]["samples"] = args.samples

        #update_config(default_config, config)
        #config = default_config

        with open(proj.new_config_file, 'w') as conf_out:
            yaml.dump(default_config, conf_out)
            print("hello")
    else:
        print("please supply a workdir")


def simulation(args):
    pass


def workflow(args):
    pass


def main():
    parser = argparse.ArgumentParser(
        prog='metapipe',
        usage='metapipe [subcommand] [options]',
        description='metapipe, a metagenomics data process pipeline'
    )
    parser.add_argument(
        '-v',
        '--version',
        action='store_true',
        default=False,
        help='print software version and exit')
    subparsers = parser.add_subparsers(
        title='available subcommands',
        metavar='')
    parser_init = subparsers.add_parser(
        'init',
        prog='metapipe init',
        description='a metapipe initialization',
        help='a metapipe initialization')
    parser_simulation = subparsers.add_parser(
        'simulation',
        prog='metapipe simulation',
        description='a simulation on metagenomics data',
        help='a simulation on metagenomics data')
    parser_workflow = subparsers.add_parser(
        'workflow',
        prog='metapipe workflow',
        description='a workflow on real metagenomics data',
        help='a workflow on real metagenomics data')

    parser_init.add_argument(
        '-q',
        '--queue',
        default='st.q',
        help='cluster queue')
    parser_init.add_argument(
        '-p',
        '--project',
        help='project id')
    parser_init.add_argument(
        '-d',
        '--workdir',
        help='project working directory')
    parser_init.add_argument(
        '-s',
        '--samples',
        help='raw fastq samples list')
    parser_init._optionals.title = 'arguments'
    parser_init.set_defaults(func=initialization)

    parser_simulation.add_argument(
        '-t',
        '--taxid',
        help='species id')
    parser_simulation.add_argument(
        '-g',
        '--genomes',
        metavar='<genomes.fasta>',
        help='genomes fasta')
    parser_simulation.add_argument(
        '-c',
        '--coverage',
        default=100,
        help='reads coverage, default: 100X')
    parser_simulation.add_argument(
        '-m',
        '--model',
        choices=['hiseq', 'novaseq', 'miseq'],
        default='hiseq',
        help='reads error model, default: hiseq')
    parser_simulation.add_argument(
        '-u',
        '--until',
        default='checkm',
        help='run step')
    parser_simulation._optionals.title = 'arguments'
    parser_simulation.set_defaults(func=simulation)

    parser_workflow.add_argument(
        '-u',
        '--until',
        default='checkm',
        help='run step')
    parser_workflow._optionals.title = 'arguments'
    parser_workflow.set_defaults(func=workflow)

    args = parser.parse_args()

    if args.version:
        print("metapipe version %s" % __version__)
        sys.exit(0)

    args.func(args)


if __name__ == '__main__':
    main()
