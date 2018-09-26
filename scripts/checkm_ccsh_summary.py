#!/usr/bin/env python
import argparse
import csv
import os
import re

"""
    script: checkm_ccs_summary.py
    exaplain:
        ccs:
            c: completeness
            c: contamination
            sh: strain heterogeneity

    a checkm completeness, contamination and
    strain heterogeneity summary for many samples
"""
__author__ = 'Jie Zhu'
__email__ = 'zhujie@genomics.cn'
__version__ = '0.1.0'
__date__ = 'Apr 16, 2018'


def get_checkm_out(checkmout_list, out_tsv):
    headers = [
        "fastq_id", "bin_number", "completeness_ge70_num",
        "completeness_ge80_num", "completeness_ge90_num",
        "completeness_ge70_rate", "completeness_ge80_rate",
        "completeness_ge90_rate", "contamination_le30_num",
        "contamination_le20_num", "contamination_le10_num",
        "contamination_le30_rate", "contamination_le20_rate",
        "contamination_le10_rate", "strain_heterogeneity_100_num",
        "strain_heterogeneity_100_rate"
    ]
    samples_bin_info = []
    with open(checkmout_list, "r") as list_handle:
        for checkmout in list_handle:
            bin_info = {}
            bin_info["fastq_id"] = os.path.basename(
                checkmout.strip()).split(".")[0]
            bin_info["bin_number"] = 0

            bin_info["completeness_ge70_num"] = 0
            bin_info["completeness_ge80_num"] = 0
            bin_info["completeness_ge90_num"] = 0

            bin_info["contamination_le30_num"] = 0
            bin_info["contamination_le20_num"] = 0
            bin_info["contamination_le10_num"] = 0

            bin_info["strain_heterogeneity_100_num"] = 0

            bin_info["completeness_ge70_rate"] = 0.00
            bin_info["completeness_ge80_rate"] = 0.00
            bin_info["completeness_ge90_rate"] = 0.00

            bin_info["contamination_le30_rate"] = 0.00
            bin_info["contamination_le20_rate"] = 0.00
            bin_info["contamination_le10_rate"] = 0.00

            bin_info["strain_heterogeneity_100_rate"] = 0.00

            with open(checkmout.strip(), 'r') as checkmout_handle:
                print("processing %s" % checkmout.strip())
                next(checkmout_handle)
                next(checkmout_handle)
                next(checkmout_handle)
                for info in checkmout_handle:
                    if not info.strip().startswith("-"):
                        info_list = re.split(r'\s+', info.strip())
                        bin_info["bin_number"] += 1

                        if float(info_list[-1]) == 100.00:
                            bin_info["strain_heterogeneity_100_num"] += 1

                        if float(info_list[-2]) <= 10.00:
                            bin_info["contamination_le10_num"] += 1
                        if float(info_list[-2]) <= 20.00:
                            bin_info["contamination_le20_num"] += 1
                        if float(info_list[-2]) <= 30.00:
                            bin_info["contamination_le30_num"] += 1

                        if float(info_list[-3]) >= 70.00:
                            bin_info["completeness_ge70_num"] += 1
                        if float(info_list[-3]) >= 80.00:
                            bin_info["completeness_ge80_num"] += 1
                        if float(info_list[-3]) >= 90.00:
                            bin_info["completeness_ge90_num"] += 1

                bin_info["completeness_ge70_rate"] = bin_info[
                    "completeness_ge70_num"] / bin_info["bin_number"]
                bin_info["completeness_ge80_rate"] = bin_info[
                    "completeness_ge80_num"] / bin_info["bin_number"]
                bin_info["completeness_ge90_rate"] = bin_info[
                    "completeness_ge90_num"] / bin_info["bin_number"]
                bin_info["contamination_le30_rate"] = bin_info[
                    "contamination_le30_num"] / bin_info["bin_number"]
                bin_info["contamination_le20_rate"] = bin_info[
                    "contamination_le20_num"] / bin_info["bin_number"]
                bin_info["contamination_le10_rate"] = bin_info[
                    "contamination_le10_num"] / bin_info["bin_number"]
                bin_info["strain_heterogeneity_100_rate"] = bin_info[
                    "strain_heterogeneity_100_num"] / bin_info["bin_number"]

            samples_bin_info.append(bin_info)

    with open(out_tsv, "w") as out_handle:
        f_csv = csv.DictWriter(out_handle, headers, delimiter="\t")
        f_csv.writeheader()
        f_csv.writerows(samples_bin_info)


def main():
    parser = argparse.ArgumentParser(description='''
        summary checkm results:
        completeness, contamination, strain heterogeneity for many samples'
        ''')
    parser.add_argument('-l', type=str, help='checkm out list')
    parser.add_argument('-o', type=str, help='output')
    args = parser.parse_args()
    get_checkm_out(args.l, args.o)


if __name__ == '__main__':
    main()
