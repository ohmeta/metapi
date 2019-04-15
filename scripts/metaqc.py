#!/usr/bin/env python

import argparse
import os
import pandas as pd

TRIM_TEMPLATE = '''fastp --in1 {raw_r1} --in2 {raw_r2} \
--out1 {trimmed_r1} --out2 {trimmed_r2} \
--compression {compression} \
{adapter_trim_params} \
--cut_front \
--cut_front_window_size {cut_front_window_size} \
--cut_front_mean_quality {cut_front_mean_quality} \
{use_slide_window} \
--n_base_limit {n_base_limit} \
--length_required {length_required} \
--thread {threads} \
--html {html} --json {json} 2> {trim_log}'''

RMHOST_BWA_TEMPLATE = '''bwa mem -t {threads} -k {seed} {host_index_base} \
{trimmed_r1} {trimmed_r2} 2> {rmhost_log} | \
tee >(samtools flagstat -@{threads} - > {flagstat}) | \
tee >(samtools stats -@{threads} - > {stat}) | \
{save_bam_params} \
samtools fastq -@{threads} -N -f 12 -F 256 -1 {rmhosted_r1} -2 {rmhosted_r2} -'''

RMHOST_BOWTIE2_TEMPLATE_BGISEQ = '''bowtie2 --threads {threads} -x {host_index_base} \
-1 {trimmed_r1} -2 {trimmed_r2} {additional_params} 2> {rmhost_log} | \
tee >(samtools flagstat -@{threads} - > {flagstat}) | \
tee >(samtools stats -@{threads} - > {stat}) | \
{save_bam_params} \
samtools view -@{threads} -SF4 - | awk -F'[/\\t]' '{{print $1}}' | sort | uniq | \
tee >(awk '{{print $0 "/1"}}' - | seqtk subseq -r {trimmed_r1} - | pigz -p {threads} -c > {rmhosted_r1}) | \
awk '{{print $0 "/2"}}' - | seqtk subseq -r {trimmed_r2} - | pigz -p {threads} -c > {rmhosted_r2}'''

RMHOST_BOWTIE2_TEMPLATE_ILLUMINA = '''bowtie2 --threads {threads} -x {host_index_base} \
-1 {trimmed_r1} -2 {trimmed_r2} {additional_params} 2> {rmhost_log} | \
tee >(samtools flagstat -@{threads} - > {flagstat}) | \
tee >(samtools stats -@{threads} - > {stat}) | \
{save_bam_params} \
samtools view -@{threads} -SF4 - | awk -F'[/\\t]' '{{print $1}}' | sort | uniq | \
tee >(seqtk subseq -r {trimmed_r1} - | pigz -p {threads} -c > {rmhosted_r1}) | \
seqtk subseq -r {trimmed_r2} - | pigz -p {threads} -c > {rmhosted_r2}'''

PLOT_BAMSTATS_TEMPLATE = '''plot-bamstats -p {prefix} {stat}'''


def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep='\s+').set_index("id", drop=False)


def get_fqpath(sample_df, sample_id, col):
    return sample_df.loc[sample_id, [col]].dropna()[0]


class trimmer:
    def __init__(self, sample_id, raw_r1, raw_r2, outdir, n_limit, min_len, adapter_trim, threads, no_use_slide_window):
        self.raw_r1 = raw_r1
        self.raw_r2 = raw_r2
        self.trimmed_r1 = os.path.join(outdir, sample_id + ".trimmed.1.fq.gz")
        self.trimmed_r2 = os.path.join(outdir, sample_id + ".trimmed.2.fq.gz")
        self.compression = 6
        self.cut_front_window_size = 4
        self.cut_front_mean_quality = 20
        if not no_use_slide_window:
            self.use_slide_window = '--cut_right --cut_right_window_size {ws} --cut_right_mean_quality {mq}'.\
                format(ws=4, mq=20)
        else:
            self.use_slide_window = '--cut_tail --cut_tail_window_size {ws} --cut_tail_mean_quality {mq}'.\
                format(ws=4, mq=20)
        self.n_base_limit = n_limit
        self.length_required = min_len
        self.threads = threads
        self.html = os.path.join(outdir, sample_id + ".fastp.html")
        self.json = os.path.join(outdir, sample_id + ".fastp.json")
        self.trim_log = os.path.join(outdir, sample_id + ".fastp.log")
        if adapter_trim:
            self.adapter_trim_params = ""
        else:
            self.adapter_trim_params = "--disable_adapter_trimming"


class rmhoster:
    def __init__(self, sample_id, seed, host_index_base, trim_dir, rmhost_dir, thread, save_bam, prefix, params=""):
        self.threads = thread
        self.seed = seed
        self.additional_params = params
        self.host_index_base = host_index_base
        self.trimmed_r1 = os.path.join(trim_dir, sample_id + ".trimmed.1.fq.gz")
        self.trimmed_r2 = os.path.join(trim_dir, sample_id + ".trimmed.2.fq.gz")
        if prefix is None:
            prefix_ = ""
        else:
            prefix_ = "." + prefix
        self.rmhosted_r1 = os.path.join(rmhost_dir, sample_id + ".rmhosted" + prefix_ + ".1.fq.gz")
        self.rmhosted_r2 = os.path.join(rmhost_dir, sample_id + ".rmhosted" + prefix_ + ".2.fq.gz")
        self.flagstat = os.path.join(rmhost_dir, sample_id + prefix_ + ".flagstat")
        self.stat = os.path.join(rmhost_dir, sample_id + prefix_ + ".alignstat")
        self.rmhost_log = os.path.join(rmhost_dir, sample_id + prefix_ + ".rmhost.log")
        if save_bam:
            self.save_bam_params = "tee >(samtools sort -@{t} -O BAM -o {sorted_bam} -) |".\
                format(t=thread, sorted_bam=os.path.join(rmhost_dir, sample_id + prefix_ + ".sorted.bam"))
        else:
            self.save_bam_params = ""


class bam_ploter:
    def __init__(self, sample_id, rmhost_dir, plot_dir, k):
        self.prefix = os.path.join(plot_dir, sample_id)
        if k is None:
            k_ = ""
        else:
            k_ = "." + k
        self.stat = os.path.join(rmhost_dir, sample_id + k_ + ".alignstat")


def main():
    parser = argparse.ArgumentParser(
        prog='metagenomics raw data quality control pipeline',
        usage='metaqc.py -s <samples.tsv> -o <output_dir> -d <host_index>',
        description='a simple pipeline to do quality control on metagenomics data')
    group = parser.add_argument_group("Options")
    group.add_argument(
        '-s',
        '--samples',
        type=str,
        help='samples.tsv')
    group.add_argument(
        '-d',
        '--database',
        type=str,
        default=None,
        help='host index base path')
    group.add_argument(
        '-m',
        '--aligner',
        type=str,
        choices=["bwa", "bowtie2"],
        help='which aligner to do rmhost')
    group.add_argument(
        '-p',
        '--platform',
        type=str,
        default="bgiseq",
        choices=['illumina', 'bgiseq'],
        help='PE reads come from which platform'
    )
    group.add_argument(
        '-t',
        '--threads',
        type=int,
        default=8,
        help='fastp and bwa, bowtie2, samtools threads')

    group = parser.add_argument_group("trimming options")
    group.add_argument(
        '-n',
        '--n_limit',
        default=5,
        type=int,
        help='max unknown base allowed'
    )
    group.add_argument(
        '-l',
        '--min_len',
        default=51,
        type=int,
        help='minimum reads length required'
    )
    group.add_argument(
        '-a',
        '--adapter_trim',
        default=False,
        action='store_true',
        help='adapter trimming, defalut: false'
    )
    group.add_argument(
        '-w',
        '--no_slide_window',
        default=False,
        action='store_true',
        help='not use slide window, default: false'
    )
    group = parser.add_argument_group("rmhost options")
    group.add_argument(
        '-b',
        '--save_bam',
        default=False,
        action='store_true',
        help='same host bam, default: false')

    group = parser.add_argument_group("bwa options(it will take effect when aligner = bwa)")
    group.add_argument(
        '-k',
        '--min_seed_len',
        default=19,
        type=int,
        help='minimum seed length for bwa mem algorithms'
    )
    group.add_argument(
        '-o',
        '--output',
        type=str,
        help='output directory')
    group.add_argument(
        '-P',
        '--prefix',
        default=None,
        type=str,
        help='output prefix'
    )
    args = parser.parse_args()

    trim_outdir = os.path.join(args.output, "00.trimmed")
    rmhost_outdir = os.path.join(args.output, "01.rmhosted")
    plot_outdir = os.path.join(args.output, "02.bam_plots")

    for i in [trim_outdir, rmhost_outdir, plot_outdir]:
        if not os.path.exists(i):
            os.makedirs(i, exist_ok=True)

    samples_df = parse_samples(args.samples)

    with open(os.path.join(args.output, "trim.sh"), 'w') as oh1, \
         open(os.path.join(args.output, "rmhost.sh"), 'w') as oh2, \
         open(os.path.join(args.output, "plotbam.sh"), 'w') as oh3:
        for sample_id in samples_df.index:
            r1 = get_fqpath(samples_df, sample_id, "fq1")
            r2 = get_fqpath(samples_df, sample_id, "fq2")

            trim_cmd = TRIM_TEMPLATE.format_map(
                vars(trimmer(sample_id, r1, r2,
                             os.path.basename(trim_outdir),
                             args.n_limit,
                             args.min_len,
                             args.adapter_trim,
                             args.threads,
                             args.no_slide_window)))
            rmhost_cmd = ""
            plotbam_cmd = ""

            if args.database is not None:
                if args.aligner == "bwa":
                    rmhost_cmd = RMHOST_BWA_TEMPLATE.format_map(
                        vars(rmhoster(sample_id,
                                      args.min_seed_len,
                                      args.database,
                                      os.path.basename(trim_outdir),
                                      os.path.basename(rmhost_outdir),
                                      args.threads,
                                      args.save_bam,
                                      args.prefix, "")))
                elif args.aligner == "bowtie2":
                    if args.platform == "bgiseq":
                        rmhost_cmd = RMHOST_BOWTIE2_TEMPLATE_BGISEQ.format_map(
                            vars(rmhoster(sample_id,
                                          args.min_seed_len,
                                          args.database,
                                          os.path.basename(trim_outdir),
                                          os.path.basename(rmhost_outdir),
                                          args.threads,
                                          args.save_bam,
                                          args.prefix, "--no-unal")))
                    elif args.platform == "illumina":
                        rmhost_cmd = RMHOST_BOWTIE2_TEMPLATE_ILLUMINA.format_map(
                            vars(rmhoster(sample_id,
                                          args.min_seed_len,
                                          args.database,
                                          os.path.basename(trim_outdir),
                                          os.path.basename(rmhost_outdir),
                                          args.threads,
                                          args.save_bam,
                                          args.prefix, "--no-unal")))
                    else:
                        print("unknown platform")

                plotbam_cmd = PLOT_BAMSTATS_TEMPLATE.format_map(
                    vars(bam_ploter(sample_id,
                                    os.path.basename(rmhost_outdir),
                                    os.path.basename(plot_outdir),
                                    args.prefix)))

            oh1.write(trim_cmd + "\n")
            oh2.write(rmhost_cmd + "\n")
            oh3.write(plotbam_cmd + '\n')

    os.chmod(os.path.join(args.output, "trim.sh"), 0o755)
    os.chmod(os.path.join(args.output, "rmhost.sh"), 0o755)
    os.chmod(os.path.join(args.output, "plotbam.sh"), 0o755)


if __name__ == '__main__':
    main()
