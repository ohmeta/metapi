#!/usr/bin/env python

import os
import sys
import glob
import subprocess


bins_dir = os.path.dirname(sys.argv[1])

bins_list_dastools = glob.glob(
    os.path.join(sys.argv[1] + "_DASTool_bins", "*." + sys.argv[2]))

if len(bins_list_dastools) > 0:
    for bin_fa in bins_list_dastools:
        binner = os.path.basename(bin_fa).split(".")[0]
        if binner != "unbinned":
            bin_fa_ = os.path.basename(bin_fa).replace(binner, binner +"_dastools")
            subprocess.run('''mv %s %s''' % (bin_fa, os.path.join(bins_dir, bin_fa_)), shell=True)