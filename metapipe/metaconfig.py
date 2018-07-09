#!/usr/bin/env python

import os
import yaml
import shutil

class config:
    '''
    project directory
    '''
    sub_dirs = ["assay", "assay/config", "results", "scripts", "sources", "study"]

    def __init__(self, work_dir):
        self.work_dir = os.path.realpath(work_dir)
        self.config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                                          "metaconfig.yaml")
    def __str__(self):
        message = "a metagenomics project has been created at {0}".format(self.work_dir)
        return message

    def create_dir(self):
        '''
        create project directory
        '''
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)
        
        for sub_dir in config.sub_dirs:
            loc = os.path.join(self.work_dir, sub_dir)
            if not os.path.exists(loc):
                os.mkdir(loc)

    def config_dir(self):
        return os.path.join(self.work_dir, "assay/config")



