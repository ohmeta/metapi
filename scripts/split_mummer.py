#!/usr/bin/env python
import argparse
import errno
import os
import shutil

from asub import submit_job
from split_fx import split_fasta


def gen_job(qry_fa, min_cluster, split_num, split_dir, job_dir, results_dir):
    nucmer = shutil.which("nucmer")
    for i in range(1, split_num + 1):
        # split/split_1.fa
        # job/mummer_1.sh
        # split/split_2.fa
        # job/mummer_2.sh
        # results/nucmer_1.delta
        job_sh = os.path.join(job_dir, "mummer_%i.sh" % (i))
        ref_fa = os.path.join(split_dir, "split_%i.fa" % (i))
        prefix = os.path.join(results_dir, "nucmer_%i" % (i))
        with open(job_sh, 'w') as job_h:
            job_h.write("%s -maxmatch -c %d %s %s -p %s\n" % (nucmer, min_cluster, ref_fa, qry_fa, prefix))

# TODO
# def merge():

def main():
    parser = argparse.ArgumentParser(description='''split reference, submit mummer array job to SGE, finally merge mummer results''')
    parser.add_argument('-ref', type=str, help='reference fasta file')
    parser.add_argument('-qry', type=str, help='query fasta file')
    parser.add_argument('-c', type=int, help='Sets the minimum length of a cluster of matches, default: 65', default=65)
    parser.add_argument('-size', type=int, help='how many seq records split into a group, default: 10000', default=10000)
    parser.add_argument('-outdir', type=str, help='output directory, default: ./', default="./")
    parser.add_argument('-queue', type=str, help='submit queue, default: st.q', default='st.q')
    parser.add_argument('-project', type=str, help='project id, default: F16ZQSB1SY2779', default='F16ZQSB1SY2779')
    parser.add_argument('-resource',type=str, help='resourse requirment, default: vf=1G,p=1', default='vf=1G,p=1')
    args = parser.parse_args()

    # make split, job, results dirs
    split_dir = os.path.join(os.path.abspath(args.outdir), "split")
    try:
        os.makedirs(split_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    job_dir = os.path.join(os.path.abspath(args.outdir), "job")
    try:
        os.makedirs(job_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    results_dir = os.path.join(os.path.abspath(args.outdir), "results")
    try:
        os.makedirs(results_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    qry_fa = os.path.abspath(args.qry)

    # split reference fasta
    split_num = split_fasta(args.ref, args.size, split_dir, True)
    gen_job(qry_fa, args.c, split_num, split_dir, job_dir, results_dir)
    submit_job("mummer", split_num, args.queue, args.project, args.resource, job_dir)

if __name__ == '__main__':
    main()
