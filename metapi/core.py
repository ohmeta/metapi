#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
from metapi import config


simulation_steps = [
    "genome_download",
    "genome_merge",
    "genome_simulate",
    "fastqc",
    "multiqc_fastqc",
    "trimming_oas1",
    "trimming_sickle",
    "trimming_fastp",
    "multiqc_fastp",
    "build_host_index_for_bwa",
    "rmhost_bwa",
    "build_host_index_for_bowtie2",
    "rmhost_bowtie2",
    "assembly_megahit",
    "assembly_idba_ud",
    "assembly_metaspades",
    "coassembly_megahit",
    "metaquast",
    "multiqc_metaquast",
    "prediction",
    "build_index_for_scaftigs",
    "align_reads_to_scaftigs",
    "coverage_metabat2",
    "coverage_maxbin2",
    "binning_metabat2",
    "binning_maxbin2",
    "filter_rename_prediction",
    "vsearch_clust_cds",
    "choose_cds_marker",
    "index_marker_cds",
    "get_marker_contigs_depth",
    "checkm_lineage_wf",
    "prokka_bins",
    "multiqc_prokka_bins",
    "metaphlan2_profilling",
    "metaphlan2_merge",
    "burst_reads",
    "all"
]

workflow_steps = [
    "sra2fq",
    "fastqc",
    "multiqc_fastqc",
    "trimming_oas1",
    "trimming_sickle",
    "trimming_fastp",
    "multiqc_fastp",
    "build_host_index_for_bwa",
    "rmhost_bwa",
    "build_host_index_for_bowtie2",
    "rmhost_bowtie2",
    "assembly_megahit",
    "assembly_idba_ud",
    "assembly_metaspades",
    "coassembly_megahit",
    "metaquast",
    "multiqc_metaquast",
    "prediction",
    "build_index_for_scaftigs",
    "align_reads_to_scaftigs",
    "coverage_metabat2",
    "coverage_maxbin2",
    "binning_metabat2",
    "binning_maxbin2",
    "filter_rename_prediction",
    "vsearch_clust_cds",
    "choose_cds_marker",
    "index_marker_cds",
    "get_marker_contigs_depth",
    "checkm_lineage_wf",
    "prokka_bins",
    "multiqc_prokka_bins",
    "metaphlan2_profilling",
    "metaphlan2_merge",
    "burst_reads",
    "all"
]


def initialization(args):
    if args.workdir:
        project = config.metaconfig(args.workdir)
        print(project.__str__())
        project.create_dirs()
        configuration, cluster = project.get_config()

        if args.begin:
            configuration["params"]["begin"] = args.begin
        if args.assembler:
            configuration["params"]["assembler"] = args.assembler
        if args.samples:
            configuration["params"]["samples"] = args.samples
        else:
            print("please supply a samples list!")
        if args.queue:
            cluster["__default__"]["queue"] = args.queue
        if args.project:
            cluster["__default__"]["project"] = args.project

        config.update_config(
            project.config_file, project.new_config_file, configuration, remove=False)
        config.update_config(
            project.cluster_file,
            project.new_cluster_file,
            cluster,
            remove=False)
    else:
        print("please supply a workdir!")
        sys.exit(1)


def simulation(args):
    if args.workdir:
        config_file = os.path.join(args.workdir, "config.yaml")
        configuration = config.parse_yaml(config_file)

        if args.taxid and not args.genomes:
            configuration["params"]["simulation"]["taxid"] = args.taxid
        if not args.taxid and args.genomes:
            configuration["params"]["simulation"]["genomes"] = args.genomes
        if args.taxid and args.genomes:
            print("can't specific taxid and genomes at same time")
        if args.n_genomes:
            configuration["params"]["simulation"]["n_genomes"] = args.n_genomes
        if args.n_reads:
            configuration["params"]["simulation"]["n_reads"] = args.n_reads
        if args.model:
            configuration["params"]["simulation"]["model"] = args.model

        config.update_config(config_file, config_file, configuration, remove=True)

        samples_df = pd.DataFrame({
            "id": ["s1", "s2", "s3"],
            "fq1": [
                os.path.join(config["results"]["raw"]["reads"], "s1_1.fq.gz"),
                os.path.join(config["results"]["raw"]["reads"], "s2_1.fq.gz"),
                os.path.join(config["results"]["raw"]["reads"], "s3_1.fq.gz")
            ],
            "fq2": [
                os.path.join(config["results"]["raw"]["reads"], "s1_2.fq.gz"),
                os.path.join(config["results"]["raw"]["reads"], "s2_2.fq.gz"),
                os.path.join(config["results"]["raw"]["reads"], "s3_2.fq.gz")
            ]
        })
        samples_df.to_csv(
            configuration["params"]["samples"],
            sep='\t',
            index=False,
            columns=["id", "fq1", "fq2"])
    else:
        print("please supply a workdir!")
        sys.exit(1)

    snakecmd = "snakemake --snakefile %s --configfile %s --until %s" % (
        configuration["snakefile"], configuration["configfile"], args.step)
    print(snakecmd)


def workflow(args):
    if args.workdir:
        config_file = os.path.join(args.workdir, "config.yaml")
        configuration = config.parse_yaml(config_file)
        if not os.path.exists(configuration["params"]["samples"]):
            print("please specific samples list on initialization step")
            sys.exit(1)
        if args.rmhost:
            configuration["params"]["rmhost"]["do"] = True
        config.update_config(config_file, config_file, configuration, remove=True)
    else:
        print("please supply a workdir!")
        sys.exit(1)

    snakecmd = "snakemake --snakefile %s --configfile %s --until %s" % (
        configuration["snakefile"], configuration["configfile"], args.step)
    print(snakecmd)


def main():
    parser = argparse.ArgumentParser(
        prog='metapi',
        usage='metapi [subcommand] [options]',
        description='metapi, a pipeline to construct a genome catalogue from metagenomics data')
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
        prog='metapi init',
        description='a metagenomics project initialization',
        help='a metagenomics project initialization')
    parser_simulation = subparsers.add_parser(
        'simulation',
        parents=[parent_parser],
        prog='metapi simulation',
        description='a simulation on metagenomics data',
        help='a simulation on metagenomics data')
    parser_workflow = subparsers.add_parser(
        'workflow',
        parents=[parent_parser],
        prog='metapi workflow',
        description='a workflow on real metagenomics data',
        help='a workflow on real metagenomics data')

    parser_init.add_argument(
        '-q', '--queue', default='st.q', help='cluster queue')
    parser_init.add_argument(
        '-p', '--project', default='st.m', help='project id')
    parser_init.add_argument('-s', '--samples', help='raw fastq samples list')
    parser_init.add_argument(
        '-b',
        '--begin',
        type=str,
        default='raw',
        choices=['raw', 'assembly'],
        help='begin to run pipeline from a specific step')
    parser_init.add_argument(
        '-a',
        '--assembler',
        nargs='*',
        metavar='<str>',
        default='metaspades',
        help='support metaspades, idba_ud, megahit')
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
        default='5M',
        help='reads coverage, default: 5M')
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
        choices=simulation_steps,
        default='checkm',
        help='run step')
    parser_simulation._optionals.title = 'arguments'
    parser_simulation.set_defaults(func=simulation)

    parser_workflow.add_argument(
        '-r',
        '--rmhost',
        action='store_true',
        default=False,
        help='need to remove host sequence? default: False')
    parser_workflow.add_argument(
        '-u',
        '--step',
        type=str,
        choices=workflow_steps,
        default='checkm',
        help='run step')
    parser_workflow._optionals.title = 'arguments'
    parser_workflow.set_defaults(func=workflow)

    args = parser.parse_args()
    try:
        if args.version:
            print("metapi version %s" % __version__)
            sys.exit(0)
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == '__main__':
    main()
