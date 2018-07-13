#!/usr/bin/env python

import os
import re
import shutil
import pandas


def parse_samples(sample_tsv):
    samples = {}
    with open(sample_tsv, 'r') as sample_h:
        next(sample_h)
        for line in sample_h:
            id, r1, r2 = re.split(r'\t+', line.strip())
            # assert r1.rstrip("_1.fq.gz") == r2.rstrip("_2.fq.gz"), "wrong paired names"
            # sample_id = os.path.basename(r1).rstrip("_1.fq.gz")
            # samples[sample_id] = [r1, r2]
            samples[id] = [r1, r2]
    return samples

def samples_df(samples_tsv):
    samples = pandas.read_table(samples_tsv).set_index("id", drop=False)
    return samples