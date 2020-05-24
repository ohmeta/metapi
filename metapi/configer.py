#!/usr/bin/env python

import argparse
import os
from ruamel.yaml import YAML


def parse_yaml(yaml_file):
    yaml = YAML()
    with open(yaml_file, "r") as f:
        return yaml.load(f)


def update_config(yaml_file_old, yaml_file_new, yaml_content, remove=True):
    yaml = YAML()
    yaml.default_flow_style = False
    if remove:
        os.remove(yaml_file_old)
    with open(yaml_file_new, "w") as f:
        yaml.dump(yaml_content, f)


class metaconfig:
    """
    config project directory
    """

    sub_dirs = [
        "assay",
        "results",
        "scripts",
        "sources",
        "study",
        "logs/00.simulate_short_reads",
        "logs/00.prepare_reads",
        "logs/00.raw_fastqc",
        "logs/00.raw_fastqc_multiqc",
        "logs/01.trimming_oas1",
        "logs/01.trimming_sickle",
        "logs/01.trimming_fastp",
        "logs/01.trimming_fastp_multiqc",
        "logs/02.rmhost_bwa_index",
        "logs/02.rmhost_bwa",
        "logs/02.rmhost_bowtie2_index",
        "logs/02.rmhost_bowtie2",
        "logs/03.qcreport",
        "logs/03.qcreport/raw",
        "logs/03.qcreport/trimming",
        "logs/03.qcreport/rmhost",
        "logs/03.qcreport_merge",
        "logs/03.qcreport_summary",
        "logs/04.assembly_megahit",
        "logs/04.assembly_idba_ud",
        "logs/04.assembly_metaspades",
        "logs/04.assembly_spades",
        "logs/04.assembly_metaquast",
        "logs/04.assembly_metaquast_multiqc",
        "logs/04.assembly_report",
        "logs/04.assembly_report_merge",
        "logs/04.coassembly_megahit",
        "logs/05.alignment_scaftigs_index",
        "logs/05.alignment_reads_scaftigs",
        "logs/05.alignment_bam_index",
        "logs/05.alignment_base_depth",
        "logs/05.alignment_report",
        "logs/06.binning_metabat2_coverage",
        "logs/06.binning_metabat2",
        "logs/06.binning_maxbin2_coverage",
        "logs/06.binning_maxbin2",
        "logs/06.binning_report",
        "logs/06.binning_report_merge",
        "logs/06.cobinning_vsearch_clust_cds",
        "logs/06.cobinning_choose_cds_marker",
        "logs/06.cobinning_index_cds_marker",
        "logs/06.cobinning_get_marker_contigs_coverage",
        "logs/07.predict_scaftigs_gene_prodigal",
        "logs/07.predict_scaftigs_gene_prokka",
        "logs/07.predict_bins_gene_prodigal",
        "logs/07.predict_bins_gene_prokka",
        "logs/08.checkm_lineage_wf",
        "logs/08.checkm_report",
        "logs/08.checkm_link_bins",
        "logs/09.classify_short_reads_kraken2",
        "logs/09.classify_hmq_bins_gtdbtk",
        "logs/10.dereplicate_drep",
        "logs/11.profiling_metaphlan2",
        "logs/11.profiling_metaphlan2_merge",
        "logs/11.profiling_jgi",
        "logs/11.profiling_jgi_merge",
        "logs/11.profiling_humann2",
        "logs/11.profiling_humann2_postprocess",
        "logs/11.profiling_humann2_join",
        "logs/11.profiling_humann2_split_straified",
        "logs/12.upload_generate_samples_info",
        "logs/12.upload_md5_short_reads",
        "logs/12.upload_generate_run_info",
        "logs/12.upload_md5_scaftigs",
        "logs/12.upload_generate_assembly_info",
    ]

    def __init__(self, work_dir):
        self.work_dir = os.path.realpath(work_dir)
        self.config_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "config.yaml"
        )
        self.cluster_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "cluster.yaml"
        )
        self.snake_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "Snakefile"
        )
        self.new_config_file = os.path.join(self.work_dir, "config.yaml")
        self.new_cluster_file = os.path.join(self.work_dir, "cluster.yaml")

    def __str__(self):
        message = """
.___  ___.  _______ .___________.    ___      .______    __
|   \/   | |   ____||           |   /   \     |   _  \  |  |
|  \  /  | |  |__   `---|  |----`  /  ^  \    |  |_)  | |  |
|  |\/|  | |   __|      |  |      /  /_\  \   |   ___/  |  |
|  |  |  | |  |____     |  |     /  _____  \  |  |      |  |
|__|  |__| |_______|    |__|    /__/     \__\ | _|      |__|

           Omics for All, Open Source for All

Thanks for using metapi.

A metagenomics project has been created at %s

Now, you can use "metapi denovo_wf":

metapi denovo_wf --list

metapi denovo_wf --run

metapi denovo_wf --debug

metapi denovo_wf --dry_run

metapi denovo_wf --qsub
""" % (
            self.work_dir
        )

        return message

    def create_dirs(self):
        """
        create project directory
        """
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)

        for sub_dir in metaconfig.sub_dirs:
            os.makedirs(os.path.join(self.work_dir, sub_dir), exist_ok=True)

    def get_config(self):
        """
        get default configuration
        """
        config = parse_yaml(self.config_file)
        cluster = parse_yaml(self.cluster_file)
        config["snakefile"] = self.snake_file
        config["configfile"] = self.new_config_file
        config["clusterfile"] = self.new_cluster_file
        return (config, cluster)


# https://github.com/Ecogenomics/CheckM/blob/master/checkm/customHelpFormatter.py
class custom_help_formatter(argparse.HelpFormatter):
    """Provide a customized format for help output.
    http://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if "%(default)" not in action.help:
            if (
                action.default != ""
                and action.default != []
                and action.default != None
                and action.default != False
            ):
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:
                        if "\n" in h:
                            lines = h.splitlines()
                            lines[0] += " (default: %(default)s)"
                            h = "\n".join(lines)
                        else:
                            h += " (default: %(default)s)"
            return h

    def _fill_text(self, text, width, indent):
        return "".join([indent + line for line in text.splitlines(True)])

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            (metavar,) = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return "%s %s" % (", ".join(parts), args_string)

            return ", ".join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()

    def _get_default_metavar_for_positional(self, action):
        return action.dest
