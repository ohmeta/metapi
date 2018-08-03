##!/usr/bin/env python
import glob
import os
import pprint
import sys

import pandas


def parse_bins(bins_dir):
    bin_list = []
    pattern = bins_dir + "/*/*bin*fa"
    for bin in glob.glob(pattern):
        bin_dict = {}
        bin_fa = os.path.basename(bin)
        bin_id = bin_fa.rstrip(".fa")
        id = ".".join(bin_fa.split(".")[:-3])
        bin_dict["bin_path"] = bin.strip()
        bin_dict["bin_id"] = bin_id
        bin_dict["id"] = id
        bin_list.append(bin_dict)
    pprint.pprint(bin_list)
    #bin_df = pandas.DataFrame(bin_list).set_index("bin_id", drop=False)
    #pprint.pprint(bin_df)
    # a = bin_df.loc["s1.bin.2", ["bin_path"]].dropna()[0]
    # print(a)


parse_bins(sys.argv[1])
