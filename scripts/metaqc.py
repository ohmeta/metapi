#!/usr/bin/env python

import argparse
import os
import pandas as pd

TRIM_TEMPLATE = '''fastp --in1 {raw_r1} --in2 {raw_r2} \
--out1 {trimmed_r1} --out2 {trimmed_r2} \
--compression {compression} \
 --disable_adapter_trimming \
 --cut_front --cut_right \
--cut_front_window_size {cut_front_window_size} \
--cut_front_mean_quality {cut_front_mean_quality} \
--cut_right_window_size {cut_right_window_size} \
--cut_right_mean_quality {cut_right_mean_quality} \
--n_base_limit {n_base_limit} \
--length_required {length_required} \
--html {html} --json {json} 2> {trim_log}'''

RMHOST_BWA_TEMPLATE = '''bwa mem -t {threads} {host_index_base} \
{trimmed_r1} {trimmed_r2} | \
tee >(samtools flagstat -@{threads} - > {flagstat}) | \
tee >(samtools stat -@{threads} - > {stat}) | \
tee >(samtools fastq -@{threads} -N -f 12 -F 256 -1 {rmhosted_r1} -2 {rmhosted_r2} -) | \
samtools sort -@%{threads} -O BAM -o {sorted_bam} - 2> {rmhost_log}'''

RMHOST_BOWTIE2_TEMPLATE = '''bowtie2 --threads {threads} -x {host_index_base} \
-1 {trimmed_r1} -2 {trimmed_r2} {additional_params} 2> {rmhost_log} | \
tee >(samtools flagstat -@{threads} - > {flagstat}) | \
tee >(samtools stat -@{threads} - > {stat}) | \
tee >(samtools sort -@{threads} -O BAM -o {sorted_bam}) | \
samtools view -@{threads} -SF4 - | awk -F'[/\\t]' '{{print $1}}' | sort | uniq | \
tee >(awk '{{print $0 "/1"}}' - | seqtk subseq -r {trimmed_r1} - | pigz -p {threads} -c > {rmhosted_r1}) | \
awk '{{print $0 "/2"}}' - | seqtk subseq -r {trimmed_r2} - | pigz -p {threads} -c > {rmhosted_r2}'''


def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep='\s+').set_index("id", drop=False)


def get_fqpath(sample_df, sample_id, col):
    return sample_df.loc[sample_id, [col]].dropna()[0]


class trimmer:
    def __init__(self, sample_id, raw_r1, raw_r2, outdir):
        self.raw_r1 = raw_r1
        self.raw_r2 = raw_r2
        self.trimmed_r1 = os.path.join(outdir, sample_id + ".trimmed.1.fq.gz")
        self.trimmed_r2 = os.path.join(outdir, sample_id + ".trimmed.2.fq.gz")
        self.compression = 6
        self.cut_front_window_size = 4
        self.cut_front_mean_quality = 20
        self.cut_right_window_size = 4
        self.cut_right_mean_quality = 20
        self.n_base_limit = 5
        self.length_required = 36
        self.html = os.path.join(outdir, sample_id + ".fastp.html")
        self.json = os.path.join(outdir, sample_id + ".fastp.json")
        self.trim_log = os.path.join(outdir, sample_id + ".fastp.log")


class rmhoster:
    def __init__(self, sample_id, host_index_base, trim_dir, rmhost_dir, params=""):
        self.threads = 8
        self.additional_params = params
        self.host_index_base = host_index_base
        self.trimmed_r1 = os.path.join(trim_dir, sample_id + ".trimmed.1.fq.gz")
        self.trimmed_r2 = os.path.join(trim_dir, sample_id + ".trimmed.2.fq.gz")
        self.rmhosted_r1 = os.path.join(rmhost_dir, sample_id + ".rmhosted.1.fq.gz")
        self.rmhosted_r2 = os.path.join(rmhost_dir, sample_id + ".rmhosted.2.fq.gz")
        self.flagstat = os.path.join(rmhost_dir, sample_id + ".flagstat")
        self.stat = os.path.join(rmhost_dir, sample_id + ".alignstat")
        self.sorted_bam = os.path.join(rmhost_dir, sample_id + ".sorted.bam")
        self.rmhost_log = os.path.join(rmhost_dir, sample_id + ".rmhost.log")


def insert_sizer():
    pass


def main():
    parser = argparse.ArgumentParser(
        prog='metagenomics raw data quality control pipeline',
        usage='metaqc.py -s <samples.tsv> -o <output_dir> -d <host_index>',
        description='a simple pipeline to calculate insert size on many samples')
    parser.add_argument(
        '-s',
        '--samples',
        type=str,
        help='samples.tsv')
    parser.add_argument(
        '-d',
        '--database',
        type=str,
        default=None,
        help='host index base path')
    parser.add_argument(
        '-m',
        '--aligner',
        type=str,
        choices=["bwa", "bowtie2"],
        help='which aligner to do rmhost')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='output directory')
    args = parser.parse_args()

    trim_outdir = os.path.join(args.output, "00.trimmed")
    rmhost_outdir = os.path.join(args.output, "01.rmhosted")

    if not os.path.exists(trim_outdir):
        os.makedirs(trim_outdir, exist_ok=True)
    if not os.path.exists(rmhost_outdir):
        os.makedirs(rmhost_outdir, exist_ok=True)

    samples_df = parse_samples(args.samples)
    for sample_id in samples_df.index:
        r1 = get_fqpath(samples_df, sample_id, "fq1")
        r2 = get_fqpath(samples_df, sample_id, "fq2")

        trim_cmd = TRIM_TEMPLATE.format_map(vars(trimmer(sample_id, r1, r2, trim_outdir)))
        rmhost_cmd = ""
        if args.database is not None:
            if args.aligner == "bwa":
                rmhost_cmd = RMHOST_BWA_TEMPLATE.format_map(
                    vars(rmhoster(sample_id, args.database, trim_outdir, rmhost_outdir)))
            elif args.aligner == "bowtie2":
                rmhost_cmd = RMHOST_BOWTIE2_TEMPLATE.format_map(
                    vars(rmhoster(sample_id, args.database, trim_outdir, rmhost_outdir, "--no-unal")))

        print(trim_cmd)
        print(rmhost_cmd)


if __name__ == '__main__':
    main()
