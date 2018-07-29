#!/usr/bin/env python

import os
import shutil

import yaml


class config:
    '''
    project directory
    '''
    sub_dirs = ["assay",
                "results", "results/00.raw",
                "scripts",
                "sources",
                "study"]

    def __init__(self, work_dir):
        self.work_dir = os.path.realpath(work_dir)
        self.config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "metaconfig.yaml")
        self.snake_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "Snakefile")
        self.new_config_file = os.path.join(self.work_dir, "metaconfig.yaml")
        self.samples_tsv = os.path.join(self.work_dir, "results/00.raw/samples.tsv")

    def __str__(self):
        message = "a metagenomics project has been created at {0}".format(
            self.work_dir)
        return message

    def create_dirs(self):
        '''
        create project directory
        '''
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)

        for sub_dir in config.sub_dirs:
            loc = os.path.join(self.work_dir, sub_dir)
            if not os.path.exists(loc):
                os.mkdir(loc)

    def get_config(self):
        '''
        get default configuration
        '''
        with open(self.config_file) as conf_in:
            conf_info = yaml.load_all(conf_in)
        return conf_info


        with open(os.path.join(self.work_dir, "metaconfig.yaml"), 'w') as conf_out:
                yaml.dump(conf_info, conf_out)
