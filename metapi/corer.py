#!/usr/bin/env python
import argparse
import os
import sys
import pandas as pd
from metapi import configer

SIMULATION_STEPS = ["genome_download", "genome_merge", "genome_simulate"]

QUALITY_CONTROL_STEPS = [
    # fastqc
    "fastqc",
    "multiqc_fastqc",
    "raw_report",
    "merge_raw_report",
    # sra2fq
    "sra2fq",
    # trimming
    "trimming_oas1",
    "trimming_sickle",
    "trimming_fastp",
    "multiqc_fastp",
    "trimming_report",
    "merge_trimming_report",
    # rmhost
    "build_host_index_for_bwa",
    "rmhost_bwa",
    "build_host_index_for_bowtie2",
    "rmhost_bowtie2",
    "rmhost_report",
    "merge_rmhost_report",
    # qc report
    "qc_report",
]

ASSEMBLY_STEPS = [
    "assembly_megahit",
    "assembly_idba_ud",
    "assembly_metaspades",
    "assembly_spades",
    "assembly_report",
    "assembly_summary",
]

COASSEMBLY_STEPS = [
    "coassembly_megahit",
    "demultiplex_kraken2_reads",
    "merge_kraken2_reads",
]

ASSEMBLY_EVULATION_STEPS = ["metaquast", "multiqc_metaquast"]

ALIGNMENT_STEPS = [
    "build_index_for_scaftigs",
    "align_reads_to_scaftigs",
    "build_index_for_bam",
    "cal_base_depth",
    "alignment_summary",
]

BINNING_STEPS = [
    "coverage_metabat2",
    "coverage_maxbin2",
    "binning_metabat2",
    "binning_maxbin2",
]

COBINNING_STEPS = [
    "filter_rename_prediction",
    "vsearch_clust_cds",
    "choose_cds_marker",
    "index_marker_cds",
    "get_marker_contigs_depth",
]

CHECKM_STEPS = [
    "checkm_lineage_wf",
    "checkm_report",
    "checkm_link_bins",
    # "checkm_coverage",
    # "checkm_profile"
]

ANNOTATION_STEPS = ["prediction", "prokka_bins", "multiqc_prokka_bins"]

CLASSIFICATION_STEPS = ["kraken2", "classification_hmq_bins_by_gtdbtk"]


DEREPLICATION_STEPS = ["drep"]

PROFILING_STEPS = [
    "metaphlan2_profiling",
    "metaphlan2_merge",
    "jgi_profiling",
    "jgi_profile_merge",
    "humann2_profiling",
    "humann2_postprocess",
    "humann2_join",
    "humann2_split_straified",
]

BURST_STEPS = [
    "burst_reads",
    # "burst_contigs",
    # "burst_bins"
]

UPLOAD_STEPS = [
    "rmhost_md5",
    "assembly_md5",
    "generate_samples_info",
    "generate_run_info",
    "generate_assembly_info",
]


SIMULATION_WORKFLOW = (
    SIMULATION_STEPS
    + QUALITY_CONTROL_STEPS
    + ASSEMBLY_STEPS
    + ASSEMBLY_EVULATION_STEPS
    + COASSEMBLY_STEPS
    + ALIGNMENT_STEPS
    + BINNING_STEPS
    + COBINNING_STEPS
    + CHECKM_STEPS
    + CLASSIFICATION_STEPS
    + DEREPLICATION_STEPS
    + ANNOTATION_STEPS
    + PROFILING_STEPS
    + BURST_STEPS
)

RUN_WORKFLOW = (
    QUALITY_CONTROL_STEPS
    + ASSEMBLY_STEPS
    + ASSEMBLY_EVULATION_STEPS
    + COASSEMBLY_STEPS
    + ALIGNMENT_STEPS
    + BINNING_STEPS
    + COBINNING_STEPS
    + CHECKM_STEPS
    + CLASSIFICATION_STEPS
    + DEREPLICATION_STEPS
    + ANNOTATION_STEPS
    + PROFILING_STEPS
    + BURST_STEPS
    + UPLOAD_STEPS
)


def init(args):
    if args.workdir:
        project = configer.metaconfig(args.workdir)
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

        configer.update_config(
            project.config_file, project.new_config_file, configuration, remove=False
        )
        configer.update_config(
            project.cluster_file, project.new_cluster_file, cluster, remove=False
        )
    else:
        print("please supply a workdir!")
        sys.exit(1)


def simulate_wf(args):
    config_file = os.path.join(args.workdir, "config.yaml")
    conf = configer.parse_yaml(config_file)
    if conf["params"]["simulate"]["do"]:
        print("Running on simulate datasets\n")
    else:
        if args.simulate:
            conf["params"]["simulate"]["do"] = True
            print("Running on simulate datasets\n")
            configer.update_config(config_file, config_file, conf, remove=True)
        else:
            print("Running on real datasets\n")

    snakecmd = "snakemake --snakefile %s --configfile %s --until %s" % (
        conf["snakefile"],
        conf["configfile"],
        args.step,
    )
    print(snakecmd)


def denovo_wf(args):
    if args.workdir:
        config_file = os.path.join(args.workdir, "config.yaml")
        configuration = configer.parse_yaml(config_file)
        if not os.path.exists(configuration["params"]["samples"]):
            print("please specific samples list on initialization step")
            sys.exit(1)
        if args.rmhost:
            configuration["params"]["rmhost"]["do"] = True
        configer.update_config(config_file, config_file, configuration, remove=True)
    else:
        print("please supply a workdir!")
        sys.exit(1)

    snakecmd = "snakemake --snakefile %s --configfile %s --until %s" % (
        configuration["snakefile"],
        configuration["configfile"],
        args.step,
    )
    print(snakecmd)


def main():
    parser = argparse.ArgumentParser(
        prog="metapi",
        usage="metapi [subcommand] [options]",
        description="metapi, a pipeline to construct a genome catalogue from metagenomics data",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-d", "--workdir", type=str, metavar="<str>", help="project workdir"
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    parser_init = subparsers.add_parser(
        "init",
        parents=[parent_parser],
        prog="metapi init",
        description="a metagenomics project initialization",
        help="a metagenomics project initialization",
    )
    parser_simulate_wf = subparsers.add_parser(
        "simulate_wf",
        parents=[parent_parser],
        prog="metapi simulation",
        description="a simulation on metagenomics data",
        help="a simulation on metagenomics data",
    )
    parser_denovo_wf = subparsers.add_parser(
        "denovo_wf",
        parents=[parent_parser],
        prog="metapi workflow",
        description="a workflow on real metagenomics data",
        help="a workflow on real metagenomics data",
    )

    parser_init.add_argument("-q", "--queue", default="st.q", help="cluster queue")
    parser_init.add_argument("-p", "--project", default="st.m", help="project id")
    parser_init.add_argument("-s", "--samples", help="raw fastq samples list")
    parser_init.add_argument(
        "-b",
        "--begin",
        type=str,
        default="raw",
        choices=["raw", "assembly"],
        help="begin to run pipeline from a specific step",
    )
    parser_init.add_argument(
        "-a",
        "--assembler",
        nargs="*",
        metavar="<str>",
        default="metaspades",
        help="support metaspades, spades, idba_ud, megahit",
    )
    parser_init._optionals.title = "arguments"
    parser_init.set_defaults(func=init)

    parser_simulate_wf.add_argument(
        "-t",
        "--taxid",
        nargs="*",
        metavar="<int>",
        help="reference database species id(sapce-separated)",
    )
    parser_simulate_wf.add_argument(
        "-g",
        "--genomes",
        type=str,
        metavar="<genomes.fasta>",
        help="genomes fasta, default: None",
    )
    parser_simulate_wf.add_argument(
        "-ng",
        "--n_genomes",
        type=int,
        metavar="<int>",
        help="genomes number, default: 6",
    )
    parser_simulate_wf.add_argument(
        "-nr",
        "--n_reads",
        type=str,
        metavar="<str>",
        default="5M",
        help="reads coverage, default: 5M",
    )
    parser_simulate_wf.add_argument(
        "-m",
        "--model",
        choices=["hiseq", "novaseq", "miseq"],
        default="hiseq",
        help="reads error model, default: hiseq",
    )
    parser_simulate_wf.add_argument(
        "-u",
        "--step",
        type=str,
        choices=SIMULATION_WORKFLOW,
        default="checkm",
        help="run step",
    )
    parser_simulate_wf._optionals.title = "arguments"
    parser_simulate_wf.set_defaults(func=simulate_wf)

    parser_denovo_wf.add_argument(
        "-r",
        "--rmhost",
        action="store_true",
        default=False,
        help="need to remove host sequence? default: False",
    )
    parser_denovo_wf.add_argument(
        "-u", "--step", type=str, choices=RUN_WORKFLOW, default="drep", help="run step",
    )
    parser_denovo_wf._optionals.title = "arguments"
    parser_denovo_wf.set_defaults(func=denovo_wf)

    args = parser.parse_args()
    try:
        if args.version:
            print("metapi version %s" % __version__)
            sys.exit(0)
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == "__main__":
    main()
