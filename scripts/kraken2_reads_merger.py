#!/usr/bin/env python

import argparse
import os
import sys
from glob import glob


def merger_reads(inputdir, outputdir):
    merger = {}
    for i in glob(inputdir.rstrip("/") + "/*/*.1.fq.gz"):
        taxid = int(os.path.basename(i).split(".")[1])
        if taxid in merger:
            merger[taxid].append(i)
        else:
            merger[taxid] = [i]

    for taxid in merger:
        r1_str = " ".join(merger[taxid])
        r2_str = r1_str.replace("1.fq.gz",  "2.fq.gz")
        r1 = os.path.join(outputdir, "%d.1.fq.gz" % taxid)
        r2 = os.path.join(outputdir, "%d.2.fq.gz" % taxid)
        if len(merger[taxid]) > 1:
            cmd = 'cat %s > %s && rm -rf %s && cat %s > %s && rm -rf %s' % (r1_str, r1, r1_str, r2_str, r2, r2_str)
            print(cmd)
        else:
            cmd = 'mv %s %s && mv %s %s' % (r1_str, r1, r2_str, r2)
            print(cmd)


def main(args_):
    parser = argparse.ArgumentParser("merge kraken2 partition reads of many samples")
    parser.add_argument(
        '-i',
        '--input_dir',
        help='a directory contains many sample-specific directory'
    )
    parser.add_argument(
        '-o',
        '--output_dir',
        help='output directory'
    )

    args = parser.parse_args(args_)
    merger_reads(args.input_dir, args.output_dir)


if __name__ == '__main__':
    main(sys.argv[1:])
