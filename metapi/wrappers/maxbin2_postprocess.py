#!/usr/bin/env python
import os
import sys
import subprocess


with os.scandir(sys.argv[1]) as itr:
    for entry in itr:
        bin_id, bin_suffix = os.path.splitext(entry.name)
        bin_name, cluster_num = bin_id.rsplit(".", maxsplit=1)
        bin_id = bin_name + "." + cluster_num.lstrip("0")
        if bin_suffix == ".fasta":
            subprocess.run('''mv %s %s''' \
                  % (os.path.join(sys.argv[1], entry.name),
                     os.path.join(sys.argv[1], bin_id + ".fa")), shell=True)


