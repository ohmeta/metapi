#!/usr/bin/env python

import glob
import os
import pprint
import sys


def run(dir_list, outdir, logdir):
    cmd_list = []
    for dir in dir_list:
        bin_list = glob.glob(dir + "/*bin*fa")
        for bin in bin_list:
            bin_id = os.path.basename(bin).rstrip(".fa")
            prokka_dir = os.path.join(outdir, bin_id)
            log = os.path.join(logdir, bin_id + ".prokka.log")
            cmd = "prokka %s --outdir %s --prefix %s --kingdom Bacteria --cpus 8 2> %s".format(
                bin.strip(), prokka_dir, bin_id, log)
            cmd_list.append(cmd)


cmd_list = run(sys.argv[1], sys.argv[2], sys.argv[3])
pprint.pprint(cmd_list)
