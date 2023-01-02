#!/usr/bin/env python

import os
import sys
import gzip
import subprocess
from Bio import SeqIO


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
                    outh.write("%s\t%f\n" %
                               (s[0], float(a) * s[1] / total_len))
                inh.close()

    args = (
        ["iss", "generate", "--cpus", str(threads), "--genomes"]
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
        subprocess.run(f"pigz -p {threads} {output_prefix}_R1.fastq", shell=True)
        subprocess.run(f"pigz -p {threads} {output_prefix}_R2.fastq", shell=True)
        os.rename(f"{output_prefix}_R1.fastq.gz", r1)
        os.rename(f"{output_prefix}_R2.fastq.gz", r2)
    else:
        sys.exit(1)


simulate_short_reads(
    input.genomes,
    snakemake.input["genomes"],
    snakemake.params["output_prefix"],
    snakemake.output["r1"],
    snakemake.output["r2"], 
    snakemake.output["abunf"],
    snakemake.params["model"],
    snakemake.params["reads_num"],
    snakemake.params["abundance"],
    snakemake.threads,
    str(snakemake.log))


