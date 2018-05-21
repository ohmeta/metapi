#!/usr/bin/env python
import argparse
import glob
import os
import re
import sys


def parse_bins_reads(bins_list, reads_list):
    bins_dict = {}
    reads_dict = {}
    with open(bins_list, 'r') as bins_h:
        for line in bins_h:
            # bin_id, bin_path = line.strip().split(" ")
            bin_id, bin_path = re.split(r'\s+|\t', line.strip())
            bins_dict[bin_id] = bin_path
    with open(reads_list, 'r') as reads_h:
        for line in reads_h:
            # reads_id, reads_1, reads_2 = line.strip().split(" ")
            reads_id, reads_1, reads_2 = re.split(r'\s+|\t', line.strip())
            reads_dict[reads_id] = (reads_1, reads_2)
    return (bins_dict, reads_dict)


def index(bins_dict, output_dir):
    index_cmd = {}
    index_dir = os.path.join(output_dir, "index")
    os.makedirs(index_dir)
    for bin_id in bins_dict:
        bin_index_dir = os.path.join(index_dir, bin_id)
        os.makedirs(bin_index_dir)
        prefix = bin_index_dir + "/" + bin_id
        cmd = "bwa index %s -p %s".format(bins_dict[bin_id], prefix)
        index_cmd[bin_id] = cmd
    return index_cmd


def mapping(reads_dict, bins_dict, output_dir):
    mapping_cmd = {}
    mapping_dir = os.path.join(output_dir, "mapping")
    os.makedirs(mapping_dir)
    for bin_id in bins_dict:
        mapping_cmd[bin_id] = []
        prefix = os.path.join(output_dir, "index/%s/%s".format(bin_id, bin_id))
        bin_mapping_dir = os.path.join(mapping_dir, bin_id)
        os.makedirs(bin_mapping_dir)
        for read_id in reads_dict:
            r1 = os.path.join(bin_mapping_dir, "%s-%s-mapped.1.fq.gz".format(
                bin_id, read_id))
            r2 = os.path.join(bin_mapping_dir, "%s-%s-mapped.2.fq.gz".format(
                bin_id, read_id))
            cmd = "bwa mem -t %s %s %s | \
                   tee >(samtools flagstat -@8 - > %s) | \
                   samtools fastq -@8 -F 12 -n -1 %s -2 %s -".format(
                prefix, reads_dict[read_id][0], reads_dict[read_id][1], r1, r2)
            mapping_cmd[bin_id].append(cmd)
    return mapping_cmd


def reassembly(assembly, bins_dict, output_dir):
    if assembly == "megahit":
        megahit(bins_dict, output_dir)
    if assembly == "spades":
        spades(bins_dict, output_dir)
    if assembly == "idba_ud":
        idba_ud(bins_dict, output_dir)
    else:
        sys.exit()


def megahit(bins_dict, output_dir):
    megahit_cmd = {}
    asm_dir = os.path.join(output_dir, "assembly")
    os.makedirs(asm_dir)
    for bin_id in bins_dict:
        bin_mapping_dir = os.path.join(output_dir, "mapping/%s".format(bin_id))
        bin_assembly_dir = os.path.join(asm_dir, bin_id)
        r1_pattern = bin_mapping_dir + "/*.1.fq.gz"
        r2_pattern = bin_mapping_dir + "/*.2.fq.gz"
        r1 = ",".join(glob.glob(r1_pattern))
        r2 = ",".join(glob.glob(r2_pattern))
        cmd = "megahit -1 %s -2 %s -t 8 --min-contig-len 500 \
               --out-dir %s --out-prefix %s --continue".format(
            r1, r2, bin_assembly_dir, bin_id)
        megahit_cmd[bin_id] = cmd
    return megahit_cmd


def spades(bins_dict, output_dir):
    print("comming soon")


def idba_ud(bins_dict, output_dir):
    print("comming soon")


def main():
    parser = argparse.ArgumentParser(description='reassembly reads')
    parser.add_argument('-rl', type=str, help='metagenomics paired reads list')
    parser.add_argument('-bl', type=str, help='bins fasta list')
    parser.add_argument(
        '-o', type=str, help='output dir, default: ./', default='./')
    parser.add_argument(
        '-index', action='store_true', help='build index, default: true')
    parser.add_argument(
        '-mapping', action='store_true', help='do mapping, default: true')
    parser.add_argument(
        '-assembly',
        type=str,
        help='standard assembly method, default: megahit',
        default='megahit')
    args = parser.parse_args()

    bins_dict, reads_dict = parse_bins_reads(args.bl, args.rl)

    if args.index:
        index_cmd = index(bins_dict, args.o)

    if args.mapping:
        mapping_cmd = mapping(reads_dict, bins_dict, args.o)

    assembly_cmd = reassembly(args.assembly, bins_dict, args.o)

    if args.index:
        with open("index.sh", 'w') as index_out, \
             open("mapping.sh", 'w') as mapping_out, \
             open("assembly.sh", 'w') as assembly_out, \
             open("all.sh", 'w') as all_out:
            for bin_id in bins_dict:
                index_out.write(index_cmd[bin_id] + "\n")
                all_out.write(index_cmd[bin_id] + "\n")
                for map_cmd in mapping_cmd[bin_id]:
                    mapping_out.write(map_cmd + "\n")
                    all_out.write(map_cmd + "\n")
                assembly_out.write(assembly_cmd[bin_id] + "\n")
                all_out.write(assembly_cmd[bin_id] + "\n")

    else:
        if args.mapping:
            with open("mapping.sh", 'w') as mapping_out, \
                 open("assembly.sh", 'w') as assembly_out, \
                 open("all.sh", 'w') as all_out:
                for bin_id in bins_dict:
                    for map_cmd in mapping_cmd[bin_id]:
                        mapping_out.write(map_cmd + "\n")
                        all_out.write(map_cmd + "\n")
                    assembly_out.write(assembly_cmd[bin_id] + "\n")
                    all_out.write(assembly_cmd[bin_id] + "\n")
        else:
            with open("assembly.sh", 'w') as assembly_out:
                for bin_id in bins_dict:
                    assembly_out.write(assembly_cmd[bin_id] + "\n")


if __name__ == "__main__":
    main()
