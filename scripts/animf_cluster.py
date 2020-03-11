#!/usr/bin/env python

import os
import pandas as pd
import drep

def check_drep_exists():
    try:
        from drep import argumentParser
        print("drep version: %s" % argumentParser.version())
    except ImportError:
        print("drep doesn't exists")


def cluster(Bdb, Cdb, work_dir):
    Ndb = pd.DataFrame()
    for bdb, name in drep.d_cluster.iteratre_clusters(Bdb, Cdb):
        genome_list = bdb["location"].tolist()
        anin_folder = os.path.join(work_dir, "ANImf_files")

        org_lengths = {}
        files = []
        deltafiles = []

        # genome1_vs_genome2.delta
        # genome1_vs_genome2.filtered.delta
        for g1 in genome_list:
            cur_folder = os.path.join(anin_folder, os.path.basename(g1))
            org_lengths[os.path.basename(g1)] = \
                drep.d_filter.calc_fasta_length(g1)
            for g2 in genome_list:
                file_name = "{0}/{1}_vs_{2}".format(
                    cur_folder,
                    os.path.basename(g1),
                    os.path.basename(g2)
                )
                deltafiles.append(file_name + ".filtered.delta")
        df = drep.d_cluster.process_deltafiles(deltafiles,
                                               org_lengths,
                                               coverage_method="larger")

        


def main():
    pass


if __name__ == '__main__':
    main()