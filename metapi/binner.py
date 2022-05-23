#!/usr/bin/env python3

import os
import subprocess
import resource
from Bio import SeqIO
import pandas as pd


def get_binning_info(bins_dir, cluster_file, assembler):
    if assembler.lower() in ["spades", "metaspades", "megahit"]:
        with os.scandir(bins_dir) as itr, open(cluster_file, "w") as oh:
            for entry in itr:
                bin_id, suffix = os.path.splitext(entry.name)
                if suffix == ".fa":
                    cluster_num = bin_id.split(".")[-1]
                    bin_fa = os.path.join(bins_dir, entry.name)
                    for seq in SeqIO.parse(bin_fa, "fasta"):
                        # graphbin
                        # oh.write("%s,%s" %
                        #         ("_".join(seq.id.split("_")[:2]), cluster_num))
                        # graphbin 2
                        oh.write(f"{seq.id},{cluster_num}\n")


def generate_bins(cluster_file, scaftigs, prefix):

    def get_accession(identifier):
        return "_".join(identifier.split("_")[:2])

    # graphbin
    # scaftigs_index = SeqIO.index(scaftigs, "fasta", key_function=get_accession)

    # graphbin2
    scaftigs_index = SeqIO.index(scaftigs, "fasta")

    df = pd.read_csv(cluster_file, names=["scaftigs_id", "bin_id"])\
           .astype({"scaftigs_id": str,
                    "bin_id": str})\
           .set_index("bin_id")

    for i in df.index.unique():
        scaftigs_id_list = df.loc[[i], "scaftigs_id"]\
                             .dropna().tolist()
        bin_fa = prefix + "." + i + ".fa"
        with open(bin_fa, 'w') as oh:
            for scaftigs_id in scaftigs_id_list:
                SeqIO.write(scaftigs_index[scaftigs_id], oh, "fasta")


def extract_bins_report(bins_report_table):
    if os.path.getsize(bins_report_table) > 0:
        bins_report = pd.read_csv(bins_report_table, sep='\t', header=[0, 1])\
                        .rename(columns={
                            "Unnamed: 0_level_1": "assembly_group",
                            "Unnamed: 1_level_1": "bin_id",
                            "Unnamed: 2_level_1": "bin_file",
                            "Unnamed: 3_level_1": "assembler",
                            "Unnamed: 4_level_1": "binner"}, level=1)

        bins_report = bins_report[[
            ("assembly_group", "assembly_group"),
            ("bin_id", "bin_id"),
            ("bin_file", "bin_file"),
            ("assembler", "assembler"),
            ("binner", "binner"),
            ("length", "sum"),
            ("length", "N50")]]
    
        bins_report.columns = ["assembly_group", "bin_id", "bin_file", "assembler", "binner", "length", "N50"]
        return bins_report
    else:
        return pd.DataFrame(columns=["assembly_group", "bin_id", "bin_file", "assembler", "binner", "length", "N50"])


'''
            table_bins = pd.read_csv(input.table_bins, sep="\t", header=[0, 1])
            table_bins = table_bins[
                [
                    ("bin_id", "Unnamed: 1_level_1"),
                    ("chr", "count"),
                    ("length", "sum"),
                    ("length", "min"),
                    ("length", "max"),
                    ("length", "std"),
                    ("length", "N50")
                ]
            ]
            table_bins.columns = [
                "user_genome",
                "contig_number",
                "contig_length_sum",
                "contig_length_min",
                "contig_length_max",
                "contig_length_std",
                "N50"
            ]
'''


def combine_jgi(jgi_list, output_file):
    #first = False
    #jgi_df_list = []
    #for jgi in input.jgi:
    #    if not first:
    #        # jgi format
    #        # contigName\tcontigLen\ttotalAvgDepth\t{sample_id}.align2combined_scaftigs.sorted.bam
    #        jgi_df_first = pd.read_csv(jgi, sep="\t")\
    #                     .loc[:, ["contigName", "contigLen", "totalAvgDepth"]]\
    #                     .dtype({"contigName": str, "contigLen": np.int32, "totalAvgDepth": np.float32})\
    #                     .set_index("contigName")
    #        jgi_df = pd.read_csv(jgi, sep="\t").iloc[:, [0, 3]]\
    #                   .dtype({"contigName": str})
    #        jgi_df[jgi_df.columns[1]] = jgi_df[jgi_df.columns[1]].astype(np.float32)
    #        jgi_df_list = [jgi_df_first, jgi_df.set_index("contigName")]
    #        first = True
    #    else:
    #        jgi_df = pd.read_csv(jgi, sep="\t").iloc[:, [0, 3]]\
    #                   .dtype({"contigName": str})
    #        jgi_df[jgi_df.columns[1]] = jgi_df[jgi_df.columns[1]].astype(np.float32)
    #        jgi_df_list.append(jgi_df.set_index("contigName"))
    ## big table, huge memory
    #pd.concat(jgi_df_list, axis=1).reset_index().to_csv(output.matrix, sep="\t", index=False)

    #matrix_list = []
    #for jgi in input.jgi:
    #    if not first:
    #        first = True
    #        with open(jgi, 'r') as ih:
    #            for line in ih:
    #                line_list = line.strip().split("\t")
    #                matrix_list.append(line_list)
    #    else:
    #        with open(jgi, 'r') as ih:
    #            count = -1
    #            for line in ih:
    #                count += 1
    #                line_list = line.strip().split("\t")
    #                matrix_list[count].append(line_list[3])

    #with open(output.matrix, 'w') as oh:
    #    for i in matrix_list:
    #        oh.write("\t".join(i) + "\n")

    # aovid OSError: Too many open files
    max_num_file = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    if len(jgi_list) > max_num_file:
        max_num_file += len(jgi_list)
        resource.setrlimit(resource.RLIMIT_NOFILE, (max_num_file, max_num_file))

    outdir = os.path.dirname(output_file)
    os.makedirs(outdir, exist_ok=True)

    files_handle = []
    for jgi in jgi_list:
        files_handle.append(open(jgi, 'r'))

    with open(output_file, 'w') as oh:
        for line in files_handle[0]:
            oh.write(line.strip())
            for handle in files_handle[1:]:
                depth = handle.readline().strip().split("\t")[3]
                oh.write(f'''\t{depth}''')
            oh.write("\n")

    for handle in files_handle:
        handle.close()