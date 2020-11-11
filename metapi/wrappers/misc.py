#!/usr/bin/env python3

import os
import argparse
import subprocess
import sys


def link_or_cat(args):
    if not os.path.exists(os.path.join(args.output_dir,
                                       args.basename + ".fq.gz")):
        if len(args.input_file) == 1:
            real_path = os.path.realpath(args.input_file[0])
            subprocess.call(
                f'''
                pushd {args.output_dir} && \
                ln -s {real_path} {args.basename}.fq.gz && \
                popd
                ''', shell=True, stdout=sys.stdout, stderr=sys.stderr)
        else:
            subprocess.call(
                f'''
                cat {args.input_file} > {args.output_dir}/{args.basename}.fq.gz
                ''', shell=True, stdout=sys.stdout, stderr=sys.stderr)


def main():
    parser = argparse.ArgumentParser("metapi misc")
    parser.add_argument("--basename", dest="basename")
    parser.add_argument("--input-file", dest="input_file", nargs="+")
    parser.add_argument("--output-dir", dest="output_dir")

    args = parser.parse_args()

    link_or_cat(args)


if __name__ == "__main__":
    main()
