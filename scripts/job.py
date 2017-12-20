#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# do something useful with the threads
threads = job_properties[threads]

# access property defined in the cluster configuration file (snakemake >=3.6.0)
job_properties["cluster"]["time"]

os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))
