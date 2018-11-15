#!/usr/bin/env python
import argparse
import os

import pandas as pd


def parse_input(flist, header=["id", "r1", "r2"]):
    df = pd.read_table(flist, sep=r'\s+', names=header, index_col=header[0])
    return df


def index(bins_df, output_dir):
    '''
    construct BWT index for each bin fasta or genome
    '''
    index_cmd_df = pd.DataFrame(columns=["cmd"])
    index_dir = os.path.join(output_dir, "00.index")
    os.makedirs(index_dir, exist_ok=True)
    for i in range(len(bins_df)):
        bin_id = bins_df.index[i]
        bin_index_dir = os.path.join(index_dir, bin_id)
        os.makedirs(bin_index_dir, exist_ok=True)
        prefix = bin_index_dir + "/" + bin_id
        cmd = "bwa index %s -p %s" % (bins_df.iloc[i]["fa"], prefix)
        index_cmd_df = pd.concat([
            index_cmd_df,
            pd.DataFrame([cmd], index=[bin_id], columns=["cmd"])
        ])
    return index_cmd_df


def mapping(reads_df, bins_df, output_dir):
    cmd_df = pd.DataFrame(columns=["cmd"])
    mapped_reads_df = pd.DataFrame(columns=["r1", "r2"])
    mapping_dir = os.path.join(output_dir, "01.mapping")
    os.makedirs(mapping_dir, exist_ok=True)
    for i in range(len(bins_df)):
        bin_id = bins_df.index[i]
        prefix = os.path.join(output_dir, "00.index/%s/%s" % (bin_id, bin_id))
        bin_mapping_dir = os.path.join(mapping_dir, bin_id)
        os.makedirs(bin_mapping_dir, exist_ok=True)
        for i in range(0, len(reads_df)):
            read_id = reads_df.index[i]
            r1 = os.path.join(bin_mapping_dir,
                              "%s-%s-mapped.r1.fq.gz" % (bin_id, read_id))
            r2 = os.path.join(bin_mapping_dir,
                              "%s-%s-mapped.r2.fq.gz" % (bin_id, read_id))
            stat = os.path.join(bin_mapping_dir,
                                "%s-%s-flagstat.txt" % (bin_id, read_id))
            cmd = "bwa mem -t 8 %s %s %s | tee >(samtools flagstat -@8 - > %s) | samtools fastq -@8 -F 12 -n -1 %s -2 %s -" % (
                prefix, reads_df.iloc[i]["r1"], reads_df.iloc[i]["r2"], stat,
                r1, r2)

            cmd_df = pd.concat([
                cmd_df,
                pd.DataFrame([[cmd]], index=[bin_id], columns=["cmd"])
            ])
            mapped_reads_df = pd.concat([
                mapped_reads_df,
                pd.DataFrame([[r1, r2]], index=[bin_id], columns=["r1", "r2"])
            ])

    return cmd_df, mapped_reads_df


def reassembly(assembler, mapped_reads_df, output_dir):
    assembly_cmd = {}
    if assembler == "megahit":
        assembly_cmd = megahit(mapped_reads_df, output_dir)
    if assembler == "spades":
        assembly_cmd = spades(mapped_reads_df, output_dir)
    if assembler == "idba_ud":
        assembly_cmd = idba_ud(mapped_reads_df, output_dir)
    return assembly_cmd


def megahit(mapped_reads_df, output_dir):
    megahit_cmd_df = pd.DataFrame(columns=["cmd"])
    asm_dir = os.path.join(output_dir, "02.assembly")
    os.makedirs(asm_dir, exist_ok=True)
    for bin_id in mapped_reads_df.index.unique():
        bin_assembly_dir = os.path.join(asm_dir, bin_id + ".megahit_out")
        r1 = ",".join(mapped_reads_df.loc[bin_id]["r1"])
        r2 = ",".join(mapped_reads_df.loc[bin_id]["r2"])
        cmd = "megahit -1 %s -2 %s -t 8 --min-contig-len 500 --out-dir %s --out-prefix %s --continue" % (
            r1, r2, bin_assembly_dir, bin_id)
        megahit_cmd_df = pd.concat([
            megahit_cmd_df,
            pd.DataFrame([[cmd]], index=[bin_id], columns=["cmd"])
        ])
    return megahit_cmd_df


def spades(mapped_dict, output_dir):
    assembly_cmd = {}
    print("comming soon")
    return assembly_cmd


def idba_ud(mapped_dict, output_dir):
    assembly_cmd = {}
    print("comming soon")
    return assembly_cmd


def main():
    parser = argparse.ArgumentParser(description='reassembly reads')
    parser.add_argument(
        '-q', '--query', type=str, help='metagenomics paired reads list')
    parser.add_argument(
        '-r', '--ref', type=str, help='bins fasta list or reference list')
    parser.add_argument(
        '-f',
        '--mapped',
        type=str,
        help='mapped reads list, together with begin=assembly')
    parser.add_argument(
        '-b',
        '--begin',
        type=str,
        default='assembly',
        choices=['alignment', 'assembly'],
        help='from a specific step to run reassembly pipeline')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default="./",
        help='output dir, default: ./')
    parser.add_argument(
        '-a',
        '--assembler',
        type=str,
        default="spades",
        help='standard assembly method, default: megahit')
    args = parser.parse_args()

    if args.begin == "alignment":
        bins_df = parse_input(args.ref, header=["id", "fa"])
        assert bins_df.index.is_unique, "bin/ref genome id is not unique"
        reads_df = parse_input(args.query, header=["id", "r1", "r2"])

        index_cmd_df = index(bins_df, args.output)
        mapping_cmd_df, mapped_reads_df = mapping(reads_df, bins_df,
                                                  args.output)
        assembly_cmd_df = reassembly(args.assembler, mapped_reads_df,
                                     args.output)

        with open("00.index.sh", 'w') as index_out, open(
                "01.mapping.sh", 'w') as mapping_out, open(
                    "02.assembly.sh", 'w') as assembly_out, open(
                        "all.sh", 'w') as all_out:
            for bin_id in bins_df.index.unique():
                index_cmd = index_cmd_df.loc[bin_id]["cmd"]
                index_out.write(index_cmd + "\n")
                all_out.write(index_cmd + "\n")
                mapping_cmd = "\n".join(mapping_cmd_df.loc[bin_id]["cmd"])
                mapping_out.write(mapping_cmd + "\n")
                all_out.write(mapping_cmd + "\n")
                assembly_cmd = assembly_cmd_df.loc[bin_id]["cmd"]
                assembly_out.write(assembly_cmd + "\n")
                all_out.write(assembly_cmd + "\n")

    elif args.begin == "assembly":
        assert args.mapped, "please supply a mapped fastq/a list"
        assert args.mapped and not (
            args.ref or args.query
        ), "don't supply ref list or query list when begin with assembly"
        mapped_reads_df = parse_input(args.mapped, header=["id", "r1", "r2"])
        assembly_cmd_df = reassembly(args.assembler, mapped_reads_df,
                                     args.output)
        with open("02.assembly.sh", 'w') as assembly_out:
            for bin_id in mapped_reads_df.index.unique():
                assembly_out.write(assembly_cmd_df.loc[bin_id]["cmd"] + "\n")


if __name__ == "__main__":
    main()
