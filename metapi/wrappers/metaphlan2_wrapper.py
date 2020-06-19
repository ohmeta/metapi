#!/usr/bin/env python

import argparse
import os
import subprocess
import sys


def run_metaphlan2(args):
    if os.path.exists(args.bowtie2out) and not os.path.exists(args.output_file):
        os.remove(args.bowtie2out)

    cmd = [
        "metaphlan2.py",
        "--tax_lev",
        args.tax_lev,
        "--min_cu_len",
        str(args.min_cu_len),
        "--stat_q",
        str(args.stat_q),
        "--stat",
        args.stat,
        "--sample_id",
        str(args.sample_id),
        "--nproc",
        str(args.nproc),
    ]
    if args.avoid_disqm:
        cmd += ["--avoid_disqm"]

    cmd_rel_ab = cmd + [
        "--index",
        args.index,
        "--bowtie2db",
        args.bowtie2db,
        "--bt2_ps",
        args.bt2_ps,
        "-t",
        "rel_ab",
        "--read_min_len",
        str(args.read_min_len),
        "--bowtie2out",
        args.bowtie2out,
        "--output_file",
        args.output_file,
        "--input_type",
        args.input_type,
        ",".join(args.input),
    ]
    cmd_rel_ab_str = " ".join(cmd_rel_ab)
    print(cmd_rel_ab_str + "\n")

    call_code = subprocess.call(
        cmd_rel_ab_str, shell=True, stdout=sys.stdout, stderr=sys.stderr
    )

    if call_code == 0:
        analysis_type = set(args.analysis_type) ^ set(["rel_ab"])
        for anat in analysis_type:
            cmd_t = cmd + [
                "-t",
                anat,
                "--output_file",
                args.output_file.replace("abundance.profile.tsv", anat + ".tsv"),
                "--input_type",
                "bowtie2out",
                args.bowtie2out,
            ]
            cmd_t_str = " ".join(cmd_t)
            print(cmd_t_str + "\n")
            call_code = subprocess.call(
                cmd_t, shell=True, stdout=sys.stdout, stderr=sys.stderr
            )
            if call_code != 0:
                print("Error: MetaPhlAn2 analysis type: %s run failed\n" % anat)
    else:
        print("Error: MetaPhlAn2 analysis type: rel_ab run failed\n")


def main():
    parser = argparse.ArgumentParser(description="MetaPhlAn2 wrapper")
    parser.add_argument(
        "--input", metavar="INPUT_FILE", type=str, nargs="+", default=None
    )
    parser.add_argument(
        "--analysis_type",
        metavar="ANALYSIS TYPE",
        type=str,
        nargs="+",
        default="rel_ab",
    )
    parser.add_argument("--input_type", type=str, choices=["multifastq", "fastq"])
    parser.add_argument("--bowtie2db", metavar="METAPHLAN_BOWTIE2_DB", type=str)
    parser.add_argument("--index", type=str, default="v20_m200")
    parser.add_argument(
        "--bt2_ps",
        metavar="BowTie2 presets",
        default="very-sensitive",
        choices=[
            "sensitive",
            "very-sensitive",
            "sensitive-local",
            "very-sensitive-local",
        ],
    )
    parser.add_argument("--bowtie2out", metavar="FILE_NAME", type=str, default=None)
    parser.add_argument(
        "--tax_lev",
        metavar="TAXONOMIC_LEVEL",
        type=str,
        choices="akpcofgst",
        default="a",
    )
    parser.add_argument("--min_cu_len", metavar="", default="2000", type=int)
    parser.add_argument("--stat_q", metavar="", type=float, default=0.1)
    parser.add_argument("--avoid_disqm", default=False, action="store_true")
    parser.add_argument(
        "--stat",
        metavar="",
        choices=["avg_g", "avg_l", "tavg_g", "tavg_l", "wavg_g", "wavg_l", "med"],
        default="tavg_g",
        type=str,
    )
    parser.add_argument("--output_file", metavar="output file", type=str)
    parser.add_argument(
        "--sample_id", metavar="value", type=str, default="Metaphlan2_Analysis"
    )
    parser.add_argument("--nproc", metavar="N", type=int, default=8)
    parser.add_argument("--read_min_len", type=int, default=70)

    args = parser.parse_args()

    run_metaphlan2(args)


if __name__ == "__main__":
    main()
