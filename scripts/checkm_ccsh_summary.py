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


def summary_init(sample_id):
    info = {}
    info["sample_id"] = sample_id
    info["bin_number"] = 0

    info["completeness_ge70_num"] = 0
    info["completeness_ge80_num"] = 0
    info["completeness_ge90_num"] = 0

    info["contamination_le30_num"] = 0
    info["contamination_le20_num"] = 0
    info["contamination_le10_num"] = 0

    info["strain_heterogeneity_le10_num"] = 0
    info["strain_heterogeneity_le20_num"] = 0
    info["strain_heterogeneity_le30_num"] = 0

    info["completeness_ge70_rate"] = 0.00
    info["completeness_ge80_rate"] = 0.00
    info["completeness_ge90_rate"] = 0.00

    info["contamination_le30_rate"] = 0.00
    info["contamination_le20_rate"] = 0.00
    info["contamination_le10_rate"] = 0.00

    info["strain_heterogeneity_le10_rate"] = 0.00
    info["strain_heterogeneity_le20_rate"] = 0.00
    info["strain_heterogeneity_le30_rate"] = 0.00

    info["high_quality_num"] = 0
    info["medium_quality_num"] = 0
    info["high_quality_rate"] = 0.00
    info["medium_quality_rate"] = 0.00

    return info


def get_checkm_out(checkmout_list, out_tsv):
    headers = [
        "sample_id", "bin_number",
        "completeness_ge90_num", "completeness_ge80_num", "completeness_ge70_num",
        "contamination_le10_num", "contamination_le20_num", "contamination_le30_num",
        "strain_heterogeneity_le10_num", "strain_heterogeneity_le20_num", "strain_heterogeneity_le30_num",
        "completeness_ge90_rate", "completeness_ge80_rate", "completeness_ge70_rate",
        "contamination_le10_rate", "contamination_le20_rate", "contamination_le30_rate",
        "strain_heterogeneity_le10_rate", "strain_heterogeneity_le20_rate", "strain_heterogeneity_le30_rate",
        "high_quality_num", "medium_quality_num",
        "high_quality_rate", "medium_quality_rate"
    ]
    samples_bin_info = []
    summary_info = summary_init("summary")

    with open(checkmout_list, "r") as list_handle:
        num = 0
        for checkmout in list_handle:
            num += 1
            bin_info = summary_init(os.path.basename(checkmout).split(".")[0])
            with open(checkmout.strip(), 'r') as checkmout_handle:
                has_bin = False
                print("processing bin %d : %s" % (num, checkmout.strip()))
                next(checkmout_handle)
                next(checkmout_handle)
                next(checkmout_handle)
                for info in checkmout_handle:
                    has_bin = True
                    if not info.strip().startswith("-"):
                        info_list = re.split(r'\s+', info.strip())
                        bin_info["bin_number"] += 1
                        summary_info["bin_number"] += 1
                        if float(info_list[-1]) <= 10.00:
                            bin_info["strain_heterogeneity_le10_num"] += 1
                            summary_info["strain_heterogeneity_le10_num"] += 1
                        if float(info_list[-1]) <= 20.00:
                            bin_info["strain_heterogeneity_le20_num"] += 1
                            summary_info["strain_heterogeneity_le20_num"] += 1
                        if float(info_list[-1]) <= 30.00:
                            bin_info["strain_heterogeneity_le30_num"] += 1
                            summary_info["strain_heterogeneity_le30_num"] += 1

                        if float(info_list[-2]) <= 10.00:
                            bin_info["contamination_le10_num"] += 1
                            summary_info["contamination_le10_num"] += 1
                        if float(info_list[-2]) <= 20.00:
                            bin_info["contamination_le20_num"] += 1
                            summary_info["contamination_le20_num"] += 1
                        if float(info_list[-2]) <= 30.00:
                            bin_info["contamination_le30_num"] += 1
                            summary_info["contamination_le30_num"] += 1

                        if float(info_list[-3]) >= 70.00:
                            bin_info["completeness_ge70_num"] += 1
                            summary_info["completeness_ge70_num"] += 1
                        if float(info_list[-3]) >= 80.00:
                            bin_info["completeness_ge80_num"] += 1
                            summary_info["completeness_ge80_num"] += 1
                        if float(info_list[-3]) >= 90.00:
                            bin_info["completeness_ge90_num"] += 1
                            summary_info["completeness_ge90_num"] += 1

                        if (float(info_list[-1]) < 0.5) and (float(info_list[-2]) < 5) and (float(info_list[-3]) > 90):
                            bin_info["high_quality_num"] += 1
                            summary_info["high_quality_num"] += 1

                        if (float(info_list[-1]) < 5) and (float(info_list[-3]) < 90) and (float(info_list[-3]) > 50):
                            bin_info["medium_quality_num"] += 1
                            summary_info["medium_quality_num"] += 1

                if has_bin:
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

                    bin_info["strain_heterogeneity_le10_rate"] = bin_info[
                                                                     "strain_heterogeneity_le10_num"] / bin_info[
                                                                     "bin_number"]
                    bin_info["strain_heterogeneity_le20_rate"] = bin_info[
                                                                     "strain_heterogeneity_le20_num"] / bin_info[
                                                                     "bin_number"]
                    bin_info["strain_heterogeneity_le30_rate"] = bin_info[
                                                                     "strain_heterogeneity_le30_num"] / bin_info[
                                                                     "bin_number"]

                    bin_info["high_quality_rate"] = bin_info["high_quality_num"] / bin_info["bin_number"]
                    bin_info["medium_quality_rate"] = bin_info["medium_quality_num"] / bin_info["bin_number"]

            samples_bin_info.append(bin_info)

    summary_info["completeness_ge70_rate"] = summary_info[
                                                 "completeness_ge70_num"] / summary_info["bin_number"]
    summary_info["completeness_ge80_rate"] = summary_info[
                                                 "completeness_ge80_num"] / summary_info["bin_number"]
    summary_info["completeness_ge90_rate"] = summary_info[
                                                 "completeness_ge90_num"] / summary_info["bin_number"]

    summary_info["contamination_le30_rate"] = summary_info[
                                                  "contamination_le30_num"] / summary_info["bin_number"]
    summary_info["contamination_le20_rate"] = summary_info[
                                                  "contamination_le20_num"] / summary_info["bin_number"]
    summary_info["contamination_le10_rate"] = summary_info[
                                                  "contamination_le10_num"] / summary_info["bin_number"]

    summary_info["strain_heterogeneity_le10_rate"] = summary_info[
                                                         "strain_heterogeneity_le10_num"] / summary_info["bin_number"]
    summary_info["strain_heterogeneity_le20_rate"] = summary_info[
                                                         "strain_heterogeneity_le20_num"] / summary_info["bin_number"]
    summary_info["strain_heterogeneity_le30_rate"] = summary_info[
                                                         "strain_heterogeneity_le30_num"] / summary_info["bin_number"]

    summary_info["high_quality_rate"] = summary_info["high_quality_num"] / summary_info["bin_number"]
    summary_info["medium_quality_rate"] = summary_info["medium_quality_num"] / summary_info["bin_number"]

    samples_bin_info.append(summary_info)

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
