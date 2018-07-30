#!/usr/bin/env python

import argparse
import os
import sys

import metaconfig

__version__ = "0.1.0"

simulation_steps = [
    "simulation", "fastqc", "trim", "rmhost", "qc_report", "assembly",
    "alignment", "binning", "checkm", "dereplication", "classification",
    "annotation"
]

workflow_steps = [
    "fastqc", "trim", "rmhost", "qc_report", "assembly", "alignment",
    "binning", "checkm", "dereplication", "classification", "annotation"
]


def initialization(args):
    if args.workdir:
        project = metaconfig.config(args.workdir)
        project.create_dirs()
        config = project.get_config()

        if args.queue:
            config["params"]["cluster"]["queue"] = args.queue
        if args.project:
            config["params"]["cluster"]["project"] = args.project
        if args.samples:
            config["results"]["raw"]["samples"] = args.samples

        metaconfig.update_config(
            project.config_file, project.new_config_file, config, remove=False)
    else:
        print("please supply a workdir!")
        sys.exit(1)


def simulation(args):
    if args.workdir:
        config_file = os.path.join(args.workdir, "metaconfig.yaml")
        config = metaconfig.parse_yaml(config_file)
        if args.taxid and not args.genomes:
            config["params"]["simulation"]["taxid"] = args.taxid
        if not args.taxid and args.genomes:
            config["params"]["simulation"]["genomes"] = args.genomes
        if args.taxid and args.genomes:
            print("can't specific taxid and genomes at same time")
        if args.n_genomes:
            config["params"]["simulation"]["n_genomes"] = args.n_genomes
        if args.n_reads:
            config["params"]["simulation"]["n_reads"] = args.n_reads
        if args.model:
            config["params"]["simulation"]["model"] = args.model
        metaconfig.update_config(config_file, config_file, config, remove=True)
    else:
        print("please supply a workdir!")
        sys.exit(1)

    snakecmd = "snakemake --snakefile %s --configfile %s --until %s" % (
        config["snakefile"], config["configfile"], args.step)
    print(snakecmd)


def workflow(args):
    pass


def main():
    parser = argparse.ArgumentParser(
        prog='metapipe',
        usage='metapipe [subcommand] [options]',
        description='metapipe, a metagenomics data process pipeline')
    parser.add_argument(
        '-v',
        '--version',
        action='store_true',
        default=False,
        help='print software version and exit')

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '-d', '--workdir', type=str, metavar='<str>', help='project workdir')

    subparsers = parser.add_subparsers(
        title='available subcommands', metavar='')
    parser_init = subparsers.add_parser(
        'init',
        parents=[parent_parser],
        prog='metapipe init',
        description='a metapipe initialization',
        help='a metapipe initialization')
    parser_simulation = subparsers.add_parser(
        'simulation',
        parents=[parent_parser],
        prog='metapipe simulation',
        description='a simulation on metagenomics data',
        help='a simulation on metagenomics data')
    parser_workflow = subparsers.add_parser(
        'workflow',
        parents=[parent_parser],
        prog='metapipe workflow',
        description='a workflow on real metagenomics data',
        help='a workflow on real metagenomics data')

    parser_init.add_argument(
        '-q', '--queue', default='st.q', help='cluster queue')
    parser_init.add_argument(
        '-p', '--project', default='st.m', help='project id')
    parser_init.add_argument('-s', '--samples', help='raw fastq samples list')
    parser_init._optionals.title = 'arguments'
    parser_init.set_defaults(func=initialization)

    parser_simulation.add_argument(
        '-t',
        '--taxid',
        nargs='*',
        metavar='<int>',
        help='reference database species id(sapce-separated)')
    parser_simulation.add_argument(
        '-g',
        '--genomes',
        type=str,
        metavar='<genomes.fasta>',
        help='genomes fasta, default: None')
    parser_simulation.add_argument(
        '-ng',
        '--n_genomes',
        type=int,
        metavar='<int>',
        help='genomes number, default: 6')
    parser_simulation.add_argument(
        '-nr',
        '--n_reads',
        type=str,
        metavar='<str>',
        default='20M',
        help='reads coverage, default: 20M')
    parser_simulation.add_argument(
        '-m',
        '--model',
        choices=['hiseq', 'novaseq', 'miseq'],
        default='hiseq',
        help='reads error model, default: hiseq')
    parser_simulation.add_argument(
        '-u',
        '--step',
        type=str,
        metavar='<str>',
        choices=simulation_steps,
        default='checkm',
        help='run step')
    parser_simulation._optionals.title = 'arguments'
    parser_simulation.set_defaults(func=simulation)

    parser_workflow.add_argument(
        '-u',
        '--step',
        type=str,
        metavar='<str>',
        default='checkm',
        help='run step')
    parser_workflow._optionals.title = 'arguments'
    parser_workflow.set_defaults(func=workflow)

    args = parser.parse_args()
    #try:
    if args.version:
        print("metapipe version %s" % __version__)
        sys.exit(0)
    args.func(args)
    #except AttributeError as e:
    #    parser.print_help()


if __name__ == '__main__':
    main()
