#!/usr/bin/env python
import argparse
import errno
import glob
import os
import re
import sys
import stat


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
    return bins_dict, reads_dict


def index(bins_dict, output_dir):
    index_cmd = {}
    index_dir = os.path.join(output_dir, "index")
    try:
        os.makedirs(index_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    for bin_id in bins_dict:
        bin_index_dir = os.path.join(index_dir, bin_id)
        try:
            os.makedirs(bin_index_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        prefix = bin_index_dir + "/" + bin_id
        cmd = "bwa index %s -p %s" % (bins_dict[bin_id], prefix)
        index_cmd[bin_id] = cmd
    return index_cmd


def mapping(reads_dict, bins_dict, output_dir):
    mapping_cmd = {}
    mapped_dict = {}
    mapping_dir = os.path.join(output_dir, "mapping")
    try:
        os.makedirs(mapping_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    for bin_id in bins_dict:
        mapping_cmd[bin_id] = []
        mapped_dict[bin_id] = {}
        mapped_dict[bin_id]["r1"] = []
        mapped_dict[bin_id]["r2"] = [] 
        prefix = os.path.join(output_dir, "index/%s/%s" % (bin_id, bin_id))
        bin_mapping_dir = os.path.join(mapping_dir, bin_id)
        try:
            os.makedirs(bin_mapping_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        for read_id in reads_dict:
            r1 = os.path.join(bin_mapping_dir,
                              "%s-%s-mapped.r1.fq.gz" % (bin_id, read_id))
            r2 = os.path.join(bin_mapping_dir,
                              "%s-%s-mapped.r2.fq.gz" % (bin_id, read_id))
            stat = os.path.join(bin_mapping_dir,
                                "%s-%s-flagstat.txt" % (bin_id, read_id))
            cmd = "bwa mem -t 8 %s %s %s | tee >(samtools flagstat -@8 - > %s) | samtools fastq -@8 -F 12 -n -1 %s -2 %s -" % (
                prefix, reads_dict[read_id][0], reads_dict[read_id][1], stat,
                r1, r2)
            mapping_cmd[bin_id].append(cmd)
            mapped_dict[bin_id]["r1"].append(r1)
            mapped_dict[bin_id]["r2"].append(r2)

    return mapping_cmd, mapped_dict


def reassembly(assembly, bins_dict, mapped_dict, output_dir):
    assembly_cmd = {}
    if assembly == "megahit":
        assembly_cmd = megahit(bins_dict, mapped_dict, output_dir)
    if assembly == "spades":
        assembly_cmd = spades(bins_dict, mapped_dict, output_dir)
    if assembly == "idba_ud":
        assembly_cmd = idba_ud(bins_dict, mapped_dict, output_dir)
    return assembly_cmd


def megahit(bins_dict, mapped_dict, output_dir):
    megahit_cmd = {}
    asm_dir = os.path.join(output_dir, "assembly")
    try:
        os.makedirs(asm_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    for bin_id in bins_dict:
        bin_assembly_dir = os.path.join(asm_dir, bin_id)
        r1 = ",".join(mapped_dict[bin_id]["r1"])
        r2 = ",".join(mapped_dict[bin_id]["r2"])
        cmd = "megahit -1 %s -2 %s -t 8 --min-contig-len 500 --out-dir %s --out-prefix %s --continue" % (
            r1, r2, bin_assembly_dir, bin_id)
        megahit_cmd[bin_id] = cmd
    return megahit_cmd


def spades(bins_dict, mapped_dict, output_dir):
    assembly_cmd = {}
    print("comming soon")
    return assembly_cmd


def idba_ud(bins_dict, mapped_dict, output_dir):
    assembly_cmd = {}
    print("comming soon")
    return assembly_cmd

def main():
    parser = argparse.ArgumentParser(description='reassembly reads')
    parser.add_argument('-rl', type=str, help='metagenomics paired reads list')
    parser.add_argument('-bl', type=str, help='bins fasta list')
    parser.add_argument(
        '-o', type=str, help='output dir, default: ./', default='./')
    parser.add_argument(
        '-assembly',
        type=str,
        help='standard assembly method, default: megahit',
        default='megahit')
    args = parser.parse_args()

    bins_dict, reads_dict = parse_bins_reads(args.bl, args.rl)
    index_cmd = index(bins_dict, args.o)
    mapping_cmd, mapped_dict = mapping(reads_dict, bins_dict, args.o)
    assembly_cmd = reassembly(args.assembly, bins_dict, mapped_dict, args.o)

    with open("index.sh", 'w') as index_out, open("mapping.sh", 'w') as mapping_out, open("assembly.sh", 'w') as assembly_out, open("all.sh", 'w') as all_out:
        for bin_id in bins_dict:
            index_out.write(index_cmd[bin_id] + "\n")
            all_out.write(index_cmd[bin_id] + "\n")
            for map_cmd in mapping_cmd[bin_id]:
                mapping_out.write(map_cmd + "\n")
                all_out.write(map_cmd + "\n")
            assembly_out.write(assembly_cmd[bin_id] + "\n")
            all_out.write(assembly_cmd[bin_id] + "\n")


if __name__ == "__main__":
    main()
