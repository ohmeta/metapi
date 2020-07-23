#!/usr/bin/env python3

import os
import argparse


def link(link_dir, batch_num, bin_list):
    bins = []
    with open(bin_list, "r") as ih:
        for line in ih:
            bins.append(os.path.abspath(line.strip()))

    os.makedirs(link_dir, exist_ok=True)

    if len(bins) > 0:
        for batch_id in range(0, len(bins), batch_num):
            batch_dir = os.path.join(link_dir, "bins_%d" % batch_id)
            os.makedirs(batch_dir, exist_ok=True)

            for bin_fa in bins[batch_id : batch_id + batch_num]:
                os.symlink(bin_fa, os.path.join(batch_dir, os.path.basename(bin_fa)))
    else:
        os.makedirs(os.path.join(link_dir, "bins_0"), exist_ok=True)


def main():
    parser = argparse.ArgumentParser("checkm link")
    parser.add_argument("--link_dir", help="a dir contains checkm input link")
    parser.add_argument(
        "--batch_num",
        type=int,
        default=500,
        help="how many bins each cehckm run, default: 500",
    )
    parser.add_argument("--bin_list", help="a file contains all bin path")
    args = parser.parse_args()

    link(args.link_dir, args.batch_num, args.bin_list)


if __name__ == "__main__":
    main()
