#!/usr/bin/env python3

import time
import gzip
from Bio import bgzf


def demultiplex(kraken2_output, r1, r2, change_seq_id, prefix, log=None):
    start_time = time.time()
    taxid_counter = {}
    demultiplexer = {}

    with open(kraken2_output, "r") as kh:
        for line in kh:
            cols = line.split("\t")
            read_id = cols[1]
            tax_name = cols[2].split("(")[0].strip()
            tax_id = int(cols[2].split("(")[-1].split(")")[0].split()[-1])

            demultiplexer[read_id] = tax_id
            if tax_id in taxid_counter:
                taxid_counter[tax_id][1] += 1
            else:
                taxid_counter[tax_id] = [tax_name, 1]
    if log is not None:
        log_h.write(
            "step_1: parse kraken2 output has spent %d s\n" % (time.time() - start_time)
        )

    start_time = time.time()
    gzip_h = {}
    for i in taxid_counter:
        gzip_h[i] = {}
        gzip_h[i]["r1"] = bgzf.BgzfWriter(prefix + ".%d.1.fq.gz" % i, "wb")
        gzip_h[i]["r2"] = bgzf.BgzfWriter(prefix + ".%d.2.fq.gz" % i, "wb")

    if r1.endswith(".gz"):
        r1_h = gzip.open(r1, "rt")
        r2_h = gzip.open(r2, "rt")
    else:
        r1_h = open(r1, "rt")
        r2_h = open(r2, "rt")

    if change_seq_id:
        sample_tag = os.path.basename(prefix)

    if log is not None:
        log_h.write("step_2: begin demultiplex taxid-reads\n")
    for read_1, read_2 in zip(r1_h, r2_h):
        read_id = read_1[1:].split("/")[0]
        if change_seq_id:
            gzip_h[demultiplexer[read_id]]["r1"].write(
                ">%s|%s%s%s%s"
                % (sample_tag, read_1[1:], next(r1_h), next(r1_h), next(r1_h))
            )
            gzip_h[demultiplexer[read_id]]["r2"].write(
                ">%s|%s%s%s%s"
                % (sample_tag, read_2[1:], next(r2_h), next(r2_h), next(r2_h))
            )
        else:
            gzip_h[demultiplexer[read_id]]["r1"].write(
                "%s%s%s%s" % (read_1, next(r1_h), next(r1_h), next(r1_h))
            )
            gzip_h[demultiplexer[read_id]]["r2"].write(
                "%s|%s%s%s%s" % (read_2, next(r2_h), next(r2_h), next(r2_h))
            )
    if log is not None:
        log_h.write(
            "step_2: demultiplex taxid-reads has spent %d s\n"
            % (time.time() - start_time)
        )
