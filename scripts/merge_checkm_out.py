#!/usr/bin/env python

import argparse
import pandas as pd
import re
from glob import glob
import sys
from pprint import pprint


def merge(checkm_list, sort_by):
    df = pd.DataFrame()
    if re.search(r'\*', checkm_list[0]):
        checkm_list_ = glob(checkm_list[0])
    else:
        checkm_list_ = checkm_list
    for checkm_file in checkm_list_:
        checkm_df = pd.DataFrame(columns=["bin_id", "marker_lineage",
                                          "genomes", "markers", "marker_sets",
                                          "0", "1", "2", "3", "4", "5+",
                                          "completeness", "contamination", "strain_heterogeneity"])
        with open(checkm_file, 'r') as ih:
            next(ih), next(ih), next(ih)
            for line in ih:
                if not line.startswith("--"):
                    line_list = re.split(r'\s+', line.strip())
                    checkm_df = checkm_df.append({"bin_id": line_list[0],
                                                  "marker_lineage": "-".join(line_list[1:3]),
                                                  "genomes": line_list[3],
                                                  "markers": line_list[4],
                                                  "marker_sets": line_list[5],
                                                  "0": line_list[6],
                                                  "1": line_list[7],
                                                  "2": line_list[8],
                                                  "3": line_list[9],
                                                  "4": line_list[10],
                                                  "5+": line_list[11],
                                                  "completeness": line_list[12],
                                                  "contamination": line_list[13],
                                                  "strain_heterogeneity": line_list[14]}, ignore_index=True)
        df = pd.concat([df, checkm_df])
        if sort_by == "completeness":
            df = df.sort_values(by=["completeness", "contamination", "strain_heterogeneity"],
                                ascending=[False, True, True])
        else:
            df = df.sort_values(by="bin_id")
    return df


def main():
    parser = argparse.ArgumentParser("merge many checkm out txt to one")
    parser.add_argument('-l', '--list', nargs='*', help='checkm out txt list, separated by spaces')
    parser.add_argument('-o', '--output', default=sys.stdout,
                        help='merge results, if not specific it, will print stdout')
    parser.add_argument('-s', '--sort', choices=['bin_id', 'completeness'], default="completeness",
                        help='sort merged checkm output')
    args = parser.parse_args()

    df = merge(args.list, args.sort)
    df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
