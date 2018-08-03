#!/usr/bin/env python

import glob
import os
import pprint
import sys


def run(dir_list, outdir, logdir):
    cmd_list = []
    count = 1
    with open(dir_list) as f:
        for dir in f:
            count += 1
            bin_list = glob.glob(dir.strip() + "/*bin*fa")
            for bin in bin_list:
                bin_id = os.path.basename(bin).rstrip(".fa")
                prokka_dir = os.path.join(outdir,
                                          os.path.basename(dir.strip()))
                log = os.path.join(logdir, bin_id + ".prokka.log")
                cmd = "prokka %s --outdir %s --prefix %s --kingdom Bacteria --cpus 8 2> %s" % (
                    bin.strip(), prokka_dir, bin_id, log)
                cmd_list.append(cmd)
            if count == 2:
                break
    return cmd_list


cmd_list = run(sys.argv[1], sys.argv[2], sys.argv[3])
pprint.pprint(cmd_list)
