#!/usr/bin/env python3

import os
import gzip
import sys
import subprocess
import pandas as pd
from Bio import SeqIO


def parse_genomes(config):
    header = ["id", "genome", "abundance", "reads_num", "model"]

    genomes_df = pd.read_csv(
        config["params"]["samples"], sep="\t"
    ).set_index("id", drop=False)

    cancel = False
    for i in header:
        if i not in genomes_df.columns:
            cancel = True
            print(
                'Error: "%s" not in %s header'
                % (i, config["params"]["simulate"]["genomes"])
            )

    for i in genomes_df.index.unique():
        if "." in i:
            cancel = True
            print('Error: sample id %s contains ".", please remove all "."' % i)

    if cancel:
        sys.exit(1)

    genomes_df["fq1"] = genomes_df.apply(
        lambda x: os.path.join(
            config["output"]["simulate"], "short_reads/%s.simulate.1.fq.gz" % x["id"],
        ),
        axis=1,
    )
    genomes_df["fq2"] = genomes_df.apply(
        lambda x: os.path.join(
            config["output"]["simulate"], "short_reads/%s.simulate.2.fq.gz" % x["id"],
        ),
        axis=1,
    )
    return genomes_df


def simulate_short_reads(
    genomes, output_prefix, r1, r2, abunf, model, reads_num, abundance, threads, logf,
):
    if len(abundance) != 0:
        with open(abunf, "w") as outh:
            for (g, a) in zip(genomes, abundance):
                inh = gzip.open(g, "rt") if g.endswith(".gz") else open(g, "r")
                genome = []
                total_len = 0
                for record in SeqIO.parse(inh, "fasta"):
                    total_len += len(record.seq)
                    genome.append((record.id, len(record.seq)))
                for s in genome:
                    outh.write("%s\t%f\n" % (s[0], float(a) * s[1] / total_len))
                inh.close()

    args = (
        ["iss", "generate", "--compress", "--cpus", str(threads), "--genomes"]
        + genomes
        + ["--n_reads", reads_num, "--model", model, "--output", output_prefix]
    )

    if len(abundance) != 0:
        args += ["--abundance_file", abunf]
    print(" ".join(args))
    env = os.environ.copy()
    proc = subprocess.Popen(
        args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env, encoding="utf-8",
    )
    output, error = proc.communicate()

    with open(logf, "w") as logh:
        logh.write(error)

    if proc.returncode == 0:
        if len(abundance) == 0:
            default_abunf = output_prefix + "_abundance.txt"
            if os.path.exists(default_abunf):
                os.rename(default_abunf, abunf)
        os.rename(output_prefix + "_R1.fastq.gz", r1)
        os.rename(output_prefix + "_R2.fastq.gz", r2)
    else:
        sys.exit(1)


def get_simulate_info(genomes_df, wildcards, col):
    return genomes_df.loc[[wildcards.sample], col].dropna().tolist()
