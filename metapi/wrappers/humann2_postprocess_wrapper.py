#!/usr/bin/env python

import argparse
import subprocess
import sys


def regroup_table(args):
    for i in range(0, len(args.groups)):
        cmd = [
            "humann2_regroup_table",
            "--input",
            args.input,
            "--groups",
            args.groups[i],
            "--function",
            args.function,
            "--output",
            args.output[i],
        ]
        cmd_str = " ".join(cmd)
        print(cmd_str + "\n")
        subprocess.call(cmd_str, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def join_tables(args):
    for i in range(0, len(args.output)):
        cmd = [
            "humann2_join_tables",
            "--input",
            args.input,
            "--output",
            args.output[i],
            "--file_name",
            args.file_name[i],
            "--search-subdirectories",
        ]
        cmd_str = " ".join(cmd)
        print(cmd_str + "\n")
        subprocess.call(cmd_str, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def split_stratified_table(args):
    for i in range(0, len(args.input)):
        cmd = [
            "humann2_split_stratified_table",
            "-i",
            args.input[i],
            "-o",
            args.output,
        ]
        cmd_str = " ".join(cmd)
        print(cmd_str + "\n")
        subprocess.call(cmd_str, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="HUMAnN2 postprocess wrapper")
    subparsers = parser.add_subparsers(title="available subcommands", metavar="")

    parser_regroup_table = subparsers.add_parser(
        "regroup_table", prog="humann2_postprocess_wrapper.py regroup_table"
    )
    parser_regroup_table.add_argument("--input", type=str)
    parser_regroup_table.add_argument("--groups", nargs="+")
    parser_regroup_table.add_argument("--output", nargs="+")
    parser_regroup_table.add_argument("--function", type=str)

    parser_regroup_table.set_defaults(func=regroup_table)

    parser_join_tables = subparsers.add_parser(
        "join_tables", prog="humann2_postprocess_wrapper.py join_tables"
    )
    parser_join_tables.add_argument("--input", type=str)
    parser_join_tables.add_argument("--output", nargs="+")
    parser_join_tables.add_argument("--file_name", nargs="+")
    parser_join_tables.set_defaults(func=join_tables)

    parser_split_stratified_table = subparsers.add_parser(
        "split_stratified_table",
        prog="humann2_postprocess_wrapper.py split_stratified_table",
    )
    parser_split_stratified_table.add_argument(
        "--input",
        nargs="+",
    )
    parser_split_stratified_table.add_argument(
        "--output",
        type=str,
    )
    parser_split_stratified_table.set_defaults(func=split_stratified_table)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
