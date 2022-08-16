#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import textwrap
from io import StringIO

import pandas as pd

import metapi

WORKFLOWS_MAG = [
    "simulate_all",
    "prepare_short_reads_all",
    "prepare_long_reads_all",
    "prepare_reads_all",
    "raw_fastqc_all",
    "raw_report_all",
    "raw_all",
    "trimming_sickle_all",
    "trimming_fastp_all",
    "trimming_report_all",
    "trimming_all",
    "rmhost_bwa_all",
    "rmhost_bowtie2_all",
    "rmhost_minimap2_all",
    "rmhost_kraken2_all",
    "rmhost_kneaddata_all",
    "rmhost_report_all",
    "rmhost_all",
    "qcreport_all",
    "assembly_megahit_all",
    "assembly_idba_ud_all",
    "assembly_metaspades_all",
    "assembly_spades_all",
    "assembly_plass_all",
    "assembly_opera_ms_all",
    "assembly_metaquast_all",
    "assembly_report_all",
    "assembly_all",
    "alignment_base_depth_all",
    "alignment_report_all",
    "alignment_all",
    "binning_metabat2_coverage_all",
    "binning_metabat2_all",
    "binning_maxbin2_all",
    "binning_concoct_all",
    "binning_graphbin2_all",
    "binning_dastools_all",
    "binning_vamb_prepare_all",
    "binning_vamb_all",
    "binning_report_all",
    "binning_all",
    "identify_virsorter2_all",
    "identify_deepvirfinder_all",
    "identify_single_all",
    "identify_phamb_all",
    "identify_multi_all",
    "identify_all",
    "predict_scaftigs_gene_prodigal_all",
    "predict_scaftigs_gene_prokka_all",
    "predict_bins_gene_prodigal_all",
    "predict_bins_gene_prokka_all",
    "predict_scaftigs_gene_all",
    "predict_bins_gene_all",
    "predict_all",
    "checkm_all",
    "dereplicate_mags_drep_all",
    "dereplicate_mags_all",
    "taxonomic_gtdbtk_all",
    "taxonomic_all",
    "upload_sequencing_all",
    "upload_assembly_all",
    "upload_all",
    "all",
]

WORKFLOWS_GENE = [
    "simulate_all",
    "prepare_reads_all",
    "raw_fastqc_all",
    "raw_report_all",
    "raw_all",
    "trimming_sickle_all",
    "trimming_fastp_all",
    "trimming_report_all",
    "trimming_all",
    "rmhost_bwa_all",
    "rmhost_bowtie2_all",
    "rmhost_minimap2_all",
    "rmhost_report_all",
    "rmhost_all",
    "qcreport_all",
    "assebmly_megahit_all",
    "assembly_idba_ud_all",
    "assembly_metaspades_all",
    "assembly_spades_all",
    "assembly_plass_all",
    "assembly_metaquast_all",
    "assembly_report_all",
    "predict_scaftigs_gene_prodigal_all",
    "predict_scaftigs_gene_prokka_all",
    "predict_scafitgs_gene_all",
    "predict_all",
    "dereplicate_gene_cdhit_all",
    "dereplicate_gene_all",
    "upload_sequencing_all",
    "upload_assembly_all",
    "upload_all",
    "all",
]


def run_snakemake(args, unknown, snakefile, workflow):
    conf = metapi.parse_yaml(args.config)

    if not os.path.exists(conf["params"]["samples"]):
        print("Please specific samples list on init step or change config.yaml manualy")
        sys.exit(1)

    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        args.config,
        "--cores",
        str(args.cores),
        "--until",
        args.task
    ] + unknown

    if "--touch" in unknown:
        pass
    elif args.conda_create_envs_only:
        cmd += ["--use-conda", "--conda-create-envs-only"]
        if args.conda_prefix is not None:
            cmd += ["--conda-prefix", args.conda_prefix]
    else:
        cmd += [
            "--rerun-incomplete",
            "--keep-going",
            "--printshellcmds",
            "--reason",
        ]

        if args.use_conda:
            cmd += ["--use-conda"]
            if args.conda_prefix is not None:
                cmd += ["--conda-prefix", args.conda_prefix]

        if args.list:
            cmd += ["--list"]
        elif args.run_local:
            cmd += ["--local-cores", str(args.local_cores),
                    "--jobs", str(args.jobs)]
        elif args.run_remote:
            cmd += ["--profile", args.profile,
                    "--local-cores", str(args.local_cores),
                    "--jobs", str(args.jobs)]
        elif args.debug:
            cmd += ["--debug-dag"]
        else: 
            cmd += ["--dry-run"]

        if args.dry_run and ("--dry-run" not in cmd):
            cmd += ["--dry-run"]
 
    cmd_str = " ".join(cmd).strip()
    print("Running metapi %s:\n%s" % (workflow, cmd_str))

    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str,
        shell=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
        env=env,
    )
    proc.communicate()

    print(f'''\nReal running cmd:\n{cmd_str}''')


def update_config_tools(conf, begin, trimmer, rmhoster, assemblers, binners):
    conf["params"]["simulate"]["do"] = False
    conf["params"]["begin"] = begin

    for trimmer_ in ["sickle", "fastp"]:
        if trimmer_ == trimmer:
            conf["params"]["trimming"][trimmer_]["do"] = True
        else:
            conf["params"]["trimming"][trimmer_]["do"] = False

    for rmhoster_ in ["bwa", "bowtie2", "minimap2", "kraken2", "kneaddata"]:
        if rmhoster_ == rmhoster:
            conf["params"]["rmhost"][rmhoster_]["do"] = True
        else:
            conf["params"]["rmhost"][rmhoster_]["do"] = False

    for assembler_ in ["megahit", "idba_ud", "metaspades", "spades", "plass", "opera_ms"]:
        if assembler_ in assemblers:
            conf["params"]["assembly"][assembler_]["do"] = True
        else:
            conf["params"]["assembly"][assembler_]["do"] = False

    for binner_ in ["metabat2", "maxbin2", "concoct", "graphbin2", "vamb", "dastools"]:
        if binner_ in binners:
            conf["params"]["binning"][binner_]["do"] = True
        else:
            conf["params"]["binning"][binner_]["do"] = False

    if begin == "simulate":
        conf["params"]["simulate"]["do"] = True
    elif begin == "rmhost":
        conf["params"]["trimming"][trimmer]["do"] = False
    elif (begin == "assembly") or (begin == "binning"):
        conf["params"]["raw"]["save_reads"] = True
        conf["params"]["raw"]["fastqc"]["do"] = False
        conf["params"]["qcreport"]["do"] = False

        conf["params"]["trimming"][trimmer]["do"] = False
        conf["params"]["rmhost"][rmhoster]["do"] = False
    return conf


def init(args, unknown):
    if args.workdir:
        project = metapi.metaconfig(args.workdir)
        print(project.__str__())
        project.create_dirs()
        conf = project.get_config()

        for env_name in conf["envs"]:
            conf["envs"][env_name] = os.path.join(os.path.realpath(args.workdir), f"envs/{env_name}.yaml")

        conf = update_config_tools(
            conf, args.begin, args.trimmer, args.rmhoster, args.assembler, args.binner
        )

        if args.samples:
            conf["params"]["samples"] = args.samples
        else:
            print("Please supply samples table")
            sys.exit(-1)

        metapi.update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )
    else:
        print("Please supply a workdir!")
        sys.exit(-1)


def mag_wf(args, unknown):
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/mag_wf.smk")
    run_snakemake(args, unknown, snakefile, "mag_wf")


def gene_wf(args, unknown):
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/gene_wf.smk")
    run_snakemake(args, unknown, snakefile, "gene_wf")


def snakemake_summary(snakefile, configfile, task):
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        configfile,
        "--until",
        task,
        "--summary",
    ]
    cmd_out = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    summary = pd.read_csv(StringIO(cmd_out.stdout.read().decode()), sep="\t")
    return summary


def sync(args, unknown):
    snakefile = os.path.join(
        os.path.dirname(__file__), f"snakefiles/{args.workflow}.smk"
    )
    summary_df = snakemake_summary(snakefile, args.config, args.task)
    conf = metapi.parse_yaml(args.config)

    if conf["params"]["simulate"]["do"]:
        samples_df = metapi.parse_genomes(
            conf["params"]["samples"], conf["output"]["simulate"], args.check_samples
        )
    else:
        samples_df = metapi.parse_samples(
            conf["params"]["samples"],
            conf["params"]["interleaved"],
            conf["params"]["reads_layout"],
            conf["params"]["begin"],
            args.check_samples,
        )

    samples_index = samples_df.index.unique()
    count = -1
    for i in range(0, len(samples_index), args.split_num):
        count += 1
        outdir = os.path.abspath(os.path.join(args.outdir, args.name + f"_{count}"))
        os.makedirs(outdir, exist_ok=True)
        samples_file = os.path.join(outdir, f"samples_{count}.tsv")
        config_file = os.path.join(outdir, "config.yaml")

        samples = samples_df.loc[
            samples_index[i : i + args.split_num],
        ]
        samples.to_csv(samples_file, sep="\t", index=False)
        conf["params"]["samples"] = samples_file
        metapi.update_config(args.config, config_file, conf, remove=False)

        summary = snakemake_summary(snakefile, config_file, args.task)
        summary_ = summary_df[summary_df.output_file.isin(summary.output_file)]
        sync_sh = os.path.join(
            args.workdir, f"sync-{args.workflow}-{args.task}_{count}.sh"
        )
        with open(sync_sh, "w") as oh:
            print(f"generating {sync_sh}")
            for output_file in summary_["output_file"]:
                if output_file != "-" and (not pd.isna(output_file)):
                    output_file_path = os.path.join(args.workdir, output_file)
                    oh.write(
                        f"rsync --archive --relative --progress {output_file_path} {outdir}\n"
                    )
            for log_file in summary_["log-file(s)"]:
                if log_file != "-" and (not pd.isna(log_file)):
                    log_file_path = os.path.join(args.workdir, log_file)
                    oh.write(
                        f"rsync --archive --relative --progress {log_file_path} {outdir}\n"
                    )

    print(f"please change current directory to {args.workdir} to run sync script")


def main():
    banner = """

  .___  ___.  _______ .___________.    ___      .______    __
  |   \/   | |   ____||           |   /   \     |   _  \  |  |
  |  \  /  | |  |__   `---|  |----`  /  ^  \    |  |_)  | |  |
  |  |\/|  | |   __|      |  |      /  /_\  \   |   ___/  |  |
  |  |  |  | |  |____     |  |     /  _____  \  |  |      |  |
  |__|  |__| |_______|    |__|    /__/     \__\ | _|      |__|

            Omics for All, Open Source for All

 A general metagenomics data mining system focus on robust microbiome research

"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(banner),
        prog="metapi",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
    )

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "-d",
        "--workdir",
        metavar="WORKDIR",
        type=str,
        default="./",
        help="project workdir",
    )
    common_parser.add_argument(
        "--check-samples",
        dest="check_samples",
        default=False,
        action="store_true",
        help="check samples, default: False",
    )

    run_parser = argparse.ArgumentParser(add_help=False)
    run_parser.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="config.yaml",
    )
    run_parser.add_argument(
        "--profile",
        type=str,
        default="./profiles/slurm",
        help="cluster profile name",
    )
    run_parser.add_argument(
        "--cores",
        type=int,
        default=240,
        help="all job cores, available on '--run-local'")
    run_parser.add_argument(
        "--local-cores",
        type=int,
        dest="local_cores",
        default=8,
        help="local job cores, available on '--run-remote'")
    run_parser.add_argument(
        "--jobs",
        type=int,
        default=30,
        help="cluster job numbers, available on '--run-remote'")
    run_parser.add_argument(
        "--list",
        default=False,
        action="store_true",
        help="list pipeline rules",
    )
    run_parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="debug pipeline",
    )
    run_parser.add_argument(
        "--dry-run",
        default=False,
        dest="dry_run",
        action="store_true",
        help="dry run pipeline",
    )
    run_parser.add_argument(
        "--run-local",
        default=False,
        action="store_true",
        help="run pipeline on local computer",
    )
    run_parser.add_argument(
        "--run-remote",
        default=False,
        action="store_true",
        help="run pipeline on remote cluster",
    )
    run_parser.add_argument(
        "--cluster-engine",
        default="slurm",
        choices=["slurm", "sge", "lsf", "pbs-torque"],
        help="cluster workflow manager engine, support slurm(sbatch) and sge(qsub)"
    )
    run_parser.add_argument("--wait", type=int, default=60, help="wait given seconds")
    run_parser.add_argument(
        "--use-conda",
        default=False,
        dest="use_conda",
        action="store_true",
        help="use conda environment",
    )
    run_parser.add_argument(
        "--conda-prefix",
        default="~/.conda/envs",
        dest="conda_prefix",
        help="conda environment prefix",
    )
    run_parser.add_argument(
        "--conda-create-envs-only",
        default=False,
        dest="conda_create_envs_only",
        action="store_true",
        help="conda create environments only",
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    parser_init = subparsers.add_parser(
        "init",
        formatter_class=metapi.custom_help_formatter,
        parents=[common_parser],
        prog="metapi init",
        help="init project",
    )
    parser_mag_wf = subparsers.add_parser(
        "mag_wf",
        formatter_class=metapi.custom_help_formatter,
        parents=[common_parser, run_parser],
        prog="metapi mag_wf",
        help="metagenome-assembly-genome pipeline",
    )
    parser_gene_wf = subparsers.add_parser(
        "gene_wf",
        formatter_class=metapi.custom_help_formatter,
        parents=[common_parser, run_parser],
        prog="metapi gene_wf",
        help="metagenome-assembly-gene pipeline",
    )

    parser_sync = subparsers.add_parser(
        "sync",
        formatter_class=metapi.custom_help_formatter,
        parents=[common_parser],
        prog="metapi sync",
        help="metapi sync project",
    )

    parser_init.add_argument(
        "-s",
        "--samples",
        type=str,
        default=None,
        help="""desired input:
samples list, tsv format required.

if begin from trimming, rmhost, or assembly:
    if it is fastq:
        the header is: [sample_id, assembly_group, binning_group, fq1, fq2]
    if it is sra:
        the header is: [sample_id, assembly_group, binning_group, sra]

if begin from simulate:
        the header is: [id, genome, abundance, reads_num, model]

""",
    )
    parser_init.add_argument(
        "-b",
        "--begin",
        type=str,
        default="trimming",
        choices=["simulate", "trimming", "rmhost", "assembly", "binning", "checkm"],
        help="pipeline starting point",
    )
    parser_init.add_argument(
        "--trimmer",
        type=str,
        default="fastp",
        required=False,
        choices=["sickle", "fastp"],
        help="which trimmer used",
    )
    parser_init.add_argument(
        "--rmhoster",
        type=str,
        default="bowtie2",
        required=False,
        choices=["bwa", "bowtie2", "minimap2", "kraken2", "kneaddata"],
        help="which rmhoster used",
    )
    parser_init.add_argument(
        "--assembler",
        nargs="+",
        default=["metaspades"],
        required=False,
        choices=["idba-ud", "megahit", "metaspades", "spades", "opera-ms"],
        help="which assembler used, required when begin with binning, can be changed in config.yaml",
    )
    parser_init.add_argument(
        "--binner",
        nargs="+",
        required=False,
        default=["metabat2", "concoct", "maxbin2", "vamb", "dastools"],
        help="wchich binner used",
    )
    parser_init.set_defaults(func=init)

    parser_mag_wf.add_argument(
        "task",
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        choices=WORKFLOWS_MAG,
        help="pipeline end point. Allowed values are " + ", ".join(WORKFLOWS_MAG),
    )
    parser_mag_wf.set_defaults(func=mag_wf)

    parser_gene_wf.add_argument(
        "task",
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        choices=WORKFLOWS_GENE,
        help="pipeline end point. Allowed values are " + ", ".join(WORKFLOWS_GENE),
    )
    parser_gene_wf.set_defaults(func=gene_wf)

    parser_sync.add_argument(
        "workflow",
        metavar="WORKFLOW",
        nargs="?",
        type=str,
        default="mag_wf",
        choices=["mag_wf", "gene_wf"],
        help="workflow. Allowed values are mag_wf, gene_wf",
    )
    parser_sync.add_argument(
        "task", 
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        help="pipeline end point",
    )
    parser_sync.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="config.yaml",
    )
    parser_sync.add_argument("--name", type=str, required=True, help="project basename")
    parser_sync.add_argument(
        "--outdir", type=str, required=True, help="sync to a directory"
    )
    parser_sync.add_argument(
        "--split-num",
        type=int,
        default=1,
        dest="split_num",
        help="split project to sync directory",
    )
    parser_sync.set_defaults(func=sync)

    args, unknown = parser.parse_known_args()

    try:
        if args.version:
            print("metapi version %s" % metapi.__version__)
            sys.exit(0)
        args.func(args, unknown)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == "__main__":
    main()
