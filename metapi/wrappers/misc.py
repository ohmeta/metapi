#!/usr/bin/env python3

import os
import argparse
import subprocess
import sys


def link_or_cat(args):
    fq_gz = os.path.join(args.output_dir, args.basename + ".fq.gz")

    if (not os.path.exists(fq_gz)) or (os.path.getsize(fq_gz) == 0):
        subprocess.call(
            f'''rm -rf {fq_gz}''',
            shell=True, stdout=sys.stdout, stderr=sys.stderr)

        if len(args.input_file) == 1:
            reads = os.path.realpath(args.input_file[0])
            subprocess.call(
                f'''
                pushd {args.output_dir} && \
                ln -s {reads} {args.basename}.fq.gz && \
                popd
                ''', shell=True, stdout=sys.stdout, stderr=sys.stderr)
        else:
            reads = " ".join(args.input_file)
            subprocess.call(
                f'''
                cat {reads} > {args.output_dir}/{args.basename}.fq.gz
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
