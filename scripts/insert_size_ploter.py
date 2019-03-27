#!/usr/bin/env python

import argparse
import pandas as pd
from glob import glob
import os
import re
from plotnine import *


def parse_bam_stats(bam_stats_list):
    insert_size_df = pd.DataFrame()
    bam_stats_list_ = []
    if re.search(r'\*', bam_stats_list[0]):
        bam_stats_list_ = glob(bam_stats_list[0])
    else:
        bam_stats_list_ = bam_stats_list

    for bam_stats_file in bam_stats_list_:
        df = pd.DataFrame(columns=["insert_size", "pairs_total",
                                   "inward_oriented_pairs",
                                   "outward_oriented_pairs",
                                   "other_pairs", "sample_id"])
        sample_id = os.path.basename(bam_stats_file).split(".")[0]
        with open(bam_stats_file, 'r') as ih:
            for line in ih:
                if line.startswith("IS"):
                    line_list = re.split(r'\s+', line.strip())
                    df = df.append({"sample_id": sample_id,
                                    "insert_size": line_list[1],
                                    "pairs_total": line_list[2],
                                    "inward_oriented_pairs": line_list[3],
                                    "outward_oriented_pairs": line_list[4],
                                    "other_pairs": line_list[5]}, ignore_index=True)

        insert_size_df = pd.concat([insert_size_df, df])
    return insert_size_df


def plot_insert_size(insert_size_df, outpdf):
    df_l = insert_size_df.melt(id_vars=["insert_size", "sample_id"],
                               value_vars=["pairs_total",
                                           "inward_oriented_pairs",
                                           "outward_oriented_pairs",
                                           "other_pairs"],
                               var_name="type",
                               value_name="count")
    is_plot = (ggplot(df_l, aes(x='insert_size', y='count'))
               + geom_point(aes(fill='type', colour='type'), size=0.2)
               + facet_wrap('~sample_id', scales='free')
               + ggtitle('insert size distribution'))
    is_plot.save(outpdf, width=16, height=16)


def main():
    parser = argparse.ArgumentParser('plot insert size for samtools bamstats')
    parser.add_argument('-i', nargs='*', help='bamstats file list, separated by spaces')
    parser.add_argument('-o', type=str, help='insert size plot output, pdf format')

    args = parser.parse_args()

    df = parse_bam_stats(args.i)
    plot_insert_size(df, args.o)


if __name__ == '__main__':
    main()



