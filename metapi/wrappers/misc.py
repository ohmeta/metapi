#!/usr/bin/env python3

import os
import argparse
import snakemake


def link_or_cat(args):
    if not os.path.exists(os.path.join(args.output_dir,
                                       args.basenmae + ".fq.gz")):
        if len(args.input) == 1:
            real_path = os.path.realpath(args.input[0])
            snakemake.shell(
                f'''
                pushd {args.output_dir} && \
                ln -s {real_path} {args.basename}.fq.gz && \
                popd
                ''')
        else:
            snakemake.shell(
                f'''
                cat {args.input} > {args.output_dir}/{args.basename}.fq.gz
                ''')


def main():
    parser = argparse.ArgumentParser("metapi misc")
    parser.add_argument("--basename")
    parser.add_argument("--input-file", dest="input_file", nargs="+")
    parser.add_argument("--output-dir", dest="output_dir")

    args = parser.parse_args()

    link_or_cat(args)


if __name__ == "__main__":
    main()
