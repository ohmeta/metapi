#!/usr/bin/env python

import argparse
import os
import shutil
import sys
import subprocess

import metaconfig
import metasample

run_steps = [
    "fastqc",
    "trim",
    "rmhost",
    "qc_report",
    "assembly",
    "alignment",
    "binning",
    "checkm",
    "dereplication",
    "classification",
    "annotation"
]
        

def snake_cmd(snakefile, configfile, step):
    snake_cmd = ""
    if step in run_steps:
        snake_cmd = "snakemake --snakefile %s --configfile %s --until %s" % (
            snakefile, configfile, step)
    else:
        print("wrong step!")
        sys.exit()
    return snake_cmd


def main():
    parser = argparse.ArgumentParser(description='metapipe')

    cluster = parser.add_argument_group("cluster", "args for sge cluster")
    cluster.add_argument('--queue', type=str, help='queue', default='st.q')
    cluster.add_argument('--project', type=str,
                         help='project id', default='nature')

    parser.add_argument('--workdir', type=str,
                        help='project work directory', default='./temp')
    parser.add_argument('--samples', type=str, help='raw fastq',
                        default="./temp/assay/00.raw/samples.tsv")

    parser.add_argument('--rmhost', action='store',
                        help='do you want to rmhost', default=True)

    parser.add_argument('--step', type=str,
                        choices=run_steps, help='run to which step')

    args = parser.parse_args()

    project = metaconfig.config(args.workdir)
    project.create_dirs()
    project.update_config(args.samples)

    snakejob = snake_cmd(project.snake_file, project.new_config_file, args.step)
    
    print(snakejob)
    #subprocess.run(snakejob)


if __name__ == '__main__':
    main()
