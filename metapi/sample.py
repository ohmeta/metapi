#!/usr/bin/env python
import glob
import os
import pandas as pd


def samples_validator(sample_df, input_type, is_pe):
    error_count = 0
    for i in sample_df.index:
        if input_type == "fastq":
            if is_pe:
                fq1 = sample_df.loc[[i], "fq1"].dropna().tolist()
                fq2 = sample_df.loc[[i], "fq2"].dropna().tolist()
                for r1, r2 in zip(fq1, fq2):
                    if (not os.path.exists(r1)) or (not os.path.exists(r2)):
                        print("error:\t%s\t%s\t%s" % (i, r1, r2))
                        error_count += 1
            else:
                fq = sample_df.loc[[i], "fq1"].dropna().tolist()
                for r in fq:
                    if not os.path.exists(r):
                        print("error:\t%s\t%s" % (i, r))
                        error_count += 1
        elif input_type == "sra":
            for sra in sample_df.loc[[i], "sra"].dropna().tolist():
                if not os.path.exists(sra):
                    print("error:\t%s\t%s" % (i, sra))
                    error_count += 1
        else:
            print("wrong input type! just support fastq or sra")
    return error_count


def parse_samples(samples_tsv, input_type, is_pe, check=True):
    samples_df = pd.read_csv(samples_tsv, sep='\s+').set_index("id", drop=False)
    if check:
        error_count = samples_validator(samples_df, input_type, is_pe)
        if error_count == 0:
            return samples_df
        else:
            print("find %d error" % error_count)
    else:
        return samples_df


def parse_bins(bins_dir):
    bin_list = []
    for bin_ in glob.glob(bins_dir + "/*/*bin*fa"):
        bin_dict = dict()
        bin_dict["path"] = bin_.strip()
        bin_dict["id"] = os.path.basename(bin_).rstrip(".fa")
        bin_list.append(bin_dict)
    bins = pd.DataFrame(bin_list).set_index("id", drop=False)
    return bins


def get_reads(sample_df, wildcards, col):
    return sample_df.loc[[wildcards.sample], col].dropna().tolist()


def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]


def get_sample_id_(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample_, [col]].dropna()[0]


def get_bin_id(bin_df, wildcards, col):
    return bin_df.loc[wildcards.bin, [col]].dropna()[0]


def parse_cobin_samples_id(query_list):
    with open(query_list, 'r') as ih:
        samples_id = [line.strip() for line in ih]
    return samples_id


def renamed_id(samples_df, wildcards):
    return samples_df.loc[[wildcards.sample], "id_2"].dropna().tolist()[0]


def demultiplex_kraken2(kraken2_output, r1, r2, change_seq_id, prefix, log=None):
    start_time = time.time()
    taxid_counter = {}
    demultiplexer = {}

    with open(kraken2_output, 'r') as kh:
        for line in kh:
            cols = line.split('\t')
            read_id = cols[1]
            tax_name = cols[2].split("(")[0].strip()
            tax_id = int(cols[2].split("(")[-1].split(")")[0].split()[-1])

            demultiplexer[read_id] = tax_id
            if tax_id in taxid_counter:
                taxid_counter[tax_id][1] += 1
            else:
                taxid_counter[tax_id] = [tax_name, 1]
    if not log is None:
        log_h.write("step_1: parse kraken2 output has spent %d s\n" % (time.time() - start_time))

    start_time = time.time()
    gzip_h = {}
    for i in taxid_counter:
        gzip_h[i] = {}
        gzip_h[i]["r1"] = bgzf.BgzfWriter(prefix + ".%d.1.fq.gz" % i, 'wb')
        gzip_h[i]["r2"] = bgzf.BgzfWriter(prefix + ".%d.2.fq.gz" % i, 'wb')

    if r1.endswith(".gz"):
        r1_h = gzip.open(r1, 'rt')
        r2_h = gzip.open(r2, 'rt')
    else:
        r1_h = open(r1, 'rt')
        r2_h = open(r2, 'rt')

    if change_seq_id:
        sample_tag = os.path.basename(prefix)

    if not log is None:
        log_h.write("step_2: begin demultiplex taxid-reads\n")
    for read_1, read_2 in zip(r1_h, r2_h):
        read_id = read_1[1:].split("/")[0]
        if change_seq_id:
            gzip_h[demultiplexer[read_id]]["r1"].write(">%s|%s%s%s%s" %
                (
                    sample_tag, read_1[1:], next(r1_h), next(r1_h), next(r1_h)
                ))
            gzip_h[demultiplexer[read_id]]["r2"].write(">%s|%s%s%s%s" %
                (
                    sample_tag, read_2[1:], next(r2_h), next(r2_h), next(r2_h)
                ))
        else:
            gzip_h[demultiplexer[read_id]]["r1"].write("%s%s%s%s" %
                (
                    read_1, next(r1_h), next(r1_h), next(r1_h)
                ))
            gzip_h[demultiplexer[read_id]]["r2"].write("%s|%s%s%s%s" %
                (
                    read_2, next(r2_h), next(r2_h), next(r2_h)
                ))
    if not log is None:
        log_h.write("step_2: demultiplex taxid-reads has spent %d s\n" % (time.time() - start_time))
