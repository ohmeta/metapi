#!/usr/bin/env python
import argparse
import csv
import os
import re


""""
c: coverage(for contigs)
p: profile(for bins)
summary for checkm coverage based many samples
summary for checkm profile based many samples
"""


def coverage_summary(clist, cout):
    """
    contigs information
    binned or unbinned
    """
    headers = [
        'fastq_id', 'contigs_num', 'contigs_len', 'contigs_binned_num',
        'contigs_binned_len', 'contigs_unbinned_num', 'contigs_unbinned_len',
        'contigs_coverage_mean', 'contigs_binned_coverage_mean',
        'contigs_unbinned_coverage_mean', 'contigs_mapped_reads_mean',
        'contigs_binned_mapped_reads_mean',
        'contigs_unbinned_mapped_reads_mean'
    ]
    coverage_info = []
    with open(clist, 'r') as list_handle:
        for cfile in list_handle:
            info = {i: 0 for i in headers[1:]}
            info['fastq_id'] = os.path.basename(cfile.strip()).split('.')[0]
            with open(cfile.strip(), 'r') as c_handle:
                print("processing %s" % cfile.strip())
                next(c_handle)
                for line in c_handle:
                    cinfo = line.strip().split('\t')
                    if cinfo[1] == 'unbinned':
                        info['contigs_unbinned_num'] += 1
                        info['contigs_unbinned_len'] += int(cinfo[2])
                        info['contigs_unbinned_coverage_mean'] += float(
                            cinfo[-2])
                        info['contigs_unbinned_mapped_reads_mean'] += int(
                            cinfo[-1])
                    else:
                        info['contigs_binned_num'] += 1
                        info['contigs_binned_len'] += int(cinfo[2])
                        info['contigs_binned_coverage_mean'] += float(
                            cinfo[-2])
                        info['contigs_binned_mapped_reads_mean'] += int(
                            cinfo[-1])
                info['contigs_num'] = info['contigs_binned_num'] + info[
                    'contigs_unbinned_num']
                info['contigs_len'] = info['contigs_binned_len'] + info[
                    'contigs_unbinned_len']
                info['contigs_coverage_mean'] = (
                    info['contigs_binned_coverage_mean'] +
                    info['contigs_unbinned_coverage_mean']
                ) / info['contigs_num']
                info['contigs_mapped_reads_mean'] = (
                    info['contigs_binned_mapped_reads_mean'] +
                    info['contigs_unbinned_mapped_reads_mean']
                ) / info['contigs_num']
                info['contigs_binned_coverage_mean'] /= info['contigs_num']
                info['contigs_unbinned_coverage_mean'] /= info['contigs_num']
                info['contigs_binned_mapped_reads_mean'] /= info['contigs_num']
                info['contigs_unbinned_mapped_reads_mean'] /= info[
                    'contigs_num']

            coverage_info.append(info)

    with open(cout, 'w') as out_handle:
        f_tsv = csv.DictWriter(out_handle, headers, delimiter='\t')
        f_tsv.writeheader()
        f_tsv.writerows(coverage_info)


def profile_summary(plist, pout):
    """
    bins information
    """
    headers = [
        'fastq_id', 'bins_num', 'binned_size(MB)', 'unbinned_size(MB)',
        'binned_size_mean', 'binned_mapped_reads', 'binned_mapped_reads_mean',
        'binned_mapped_rate', 'binned_mapped_rate_mean',
        'unbinned_mapped_reads', 'unbinned_mapped_rate',
        'binned_populations_mean', 'binned_community_mean',
        'unbinned_community'
    ]
    profile_info = []
    with open(plist, 'r') as list_handle:
        for pfile in list_handle:
            info = {i: 0 for i in headers[1:]}
            info['fastq_id'] = os.path.basename(pfile.strip()).split('.')[0]
            with open(pfile.strip(), 'r') as p_handle:
                print("processing %s" % pfile.strip())
                next(p_handle)
                next(p_handle)
                next(p_handle)
                for line in p_handle:
                    pinfo = re.split(r'\s+', line.strip())
                    if not pinfo[0].startswith('-') and pinfo[0] != 'unbinned':
                        info['bins_num'] += 1
                        info['binned_size(MB)'] += float(pinfo[1])
                        info['binned_mapped_reads'] += int(pinfo[2])
                        info['binned_mapped_rate'] += float(pinfo[3])
                        info['binned_populations_mean'] += float(pinfo[-2])
                        info['binned_community_mean'] += float(pinfo[-1])
                    if pinfo[0] == 'unbinned':
                        info['unbinned_size(MB)'] = float(pinfo[1])
                        info['unbinned_mapped_reads'] += int(pinfo[2])
                        info['unbinned_mapped_rate'] += float(pinfo[3])
                        info['unbinned_community'] += float(pinfo[-1])
                info['binned_size_mean'] = info['binned_size(MB)'] / info[
                    'bins_num']
                info['binned_mapped_rate_mean'] = info[
                    'binned_mapped_rate'] / info['bins_num']
                info['binned_mapped_reads_mean'] = info[
                    'binned_mapped_reads'] / info['bins_num']
            profile_info.append(info)

    with open(pout, 'w') as out_handle:
        f_tsv = csv.DictWriter(out_handle, headers, delimiter='\t')
        f_tsv.writeheader()
        f_tsv.writerows(profile_info)


def main():
    parser = argparse.ArgumentParser(
        description='summary for checkm coverage and profile based many samples'
    )
    parser.add_argument(
        '-clist',
        type=str,
        help='a file contain many samples coverage file path')
    parser.add_argument(
        '-plist',
        type=str,
        help='a file contain many samples profile file path')
    parser.add_argument(
        '-cout',
        type=str,
        help='summary coverage for checkm coverage based many samples')
    parser.add_argument(
        '-pout',
        type=str,
        help='summary profile for checkm profile based many samples')
    args = parser.parse_args()

    coverage_summary(args.clist, args.cout)
    profile_summary(args.plist, args.pout)


if __name__ == '__main__':
    main()
