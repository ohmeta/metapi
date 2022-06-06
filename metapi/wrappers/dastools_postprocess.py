#!/usr/bin/env python

import os
import sys
import glob
import subprocess


bins_prefix = sys.argv[1].replace("dastools.bin", "")

bins_list = glob.glob(os.path.join(sys.argv[1] + "_DASTool_bins", "*.fa"))

if len(bins_list) > 0:
    for bin_fa in bins_list:
        if (os.path.getsize(bin_fa) > 0) and  (not "*" in bin_fa):
            binner = os.path.basename(bin_fa).split(".")[0]
            if (binner != "unbinned") and (binner != "*"):
                bin_fa_ = bins_prefix + os.path.basename(bin_fa).replace(binner, binner +"_dastools")
                subprocess.run(f'''mv {bin_fa} {bin_fa_}''', shell=True)