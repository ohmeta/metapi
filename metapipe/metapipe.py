#!/usr/bin/env python

import argparse
import metaconfig
import shutil
import yaml

data_process_protocol = [
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

def config_init(work_dir):
    project = metaconfig.config(work_dir)
    project.create_dir()
    shutil.copy2(project.config_file, project.config_dir())

'''
def write_snakefile():
'''

def main():
    parser = argparse.ArgumentParser(description='metapipe')
    parser.add_argument('--queue', type=str, help='queue', default='st.q')
    parser.add_argument('--project', type=str, help='project id', default='nature')
    parser.add_argument('--workdir', type=str, help='project work directory', default='./temp')
    parser.add_argument('--rmhost', action='store', help='do you want to rmhost', default=True)
    parser.add_argument('--run', type=str, help='local or sge cluster', default='sge')
    parser.add_argument('--step', type=str, choices=data_process_protocol, help='run to which step')
    args = parser.parse_args()

    config_init(args.workdir)

if __name__ == '__main__':
    main()