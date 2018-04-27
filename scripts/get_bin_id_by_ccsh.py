#!/usr/bin/env python
import argparse
import csv
import os
import re


def get_bin_id(checkmout_list, out_tsv, completeness, contamination):
    headers = [
        "sample_id", "bin_id", "marker_lineage", "genomes", "markers",
        "marker_sets", "completeness", "contamination", "strain_heterogeneity"
    ]
    samples_bin_info = []
    with open(checkmout_list, "r") as list_handle:
        for checkmout in list_handle:
            with open(checkmout.strip(), 'r') as checkmout_handle:
                print("processing %s" % checkmout.strip())
                sample_id = os.path.basename(checkmout.strip()).split('.')[0]
                next(checkmout_handle)
                next(checkmout_handle)
                next(checkmout_handle)
                for info in checkmout_handle:
                    if info.strip().startswith("R0"):
                        info_l = re.split(r'\s+', info.strip())
                        if (float(info_l[-2]) < contamination) and (float(
                                info_l[-3]) > completeness):
                            bin_info = {}
                            bin_info['sample_id'] = sample_id
                            bin_info["bin_id"] = info_l[0]
                            bin_info[
                                "marker_lineage"] = info_l[1] + " " + info_l[2]
                            bin_info["genomes"] = info_l[3]
                            bin_info["markers"] = info_l[4]
                            bin_info["marker_sets"] = info_l[5]
                            bin_info["completeness"] = info_l[-3]
                            bin_info["contamination"] = info_l[-2]
                            bin_info["strain_heterogeneity"] = info_l[-1]
                            samples_bin_info.append(bin_info)
    with open(out_tsv, 'w') as out_handle:
        f_tsv = csv.DictWriter(out_handle, headers, delimiter="\t")
        f_tsv.writeheader()
        f_tsv.writerows(samples_bin_info)


def main():
    parser = argparse.ArgumentParser(
        description='''get bin id by completeness cutoff and contamination
        cutoff''')
    parser.add_argument('-l', type=str, help='checkmout list of many samples')
    parser.add_argument(
        '-o',
        type=str,
        help='bin id and completeness, contamination output file')
    parser.add_argument(
        '-c1', type=float, help='completeness cutoff', default=70.0)
    parser.add_argument(
        '-c2', type=float, help='contamination cutoff', default=30.0)
    args = parser.parse_args()
    get_bin_id(args.l, args.o, args.c1, args.c2)


if __name__ == '__main__':
    main()
