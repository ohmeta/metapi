#!/usr/bin/env
import shutil
import os
from datetime import datetime
import argparse

from asub import submit_job

# TODO
#def assembly(fqlist):

def coassembly(fqlist, thread, outdir, prefix, queue, project, resource):
    r1 = []
    r2 = []
    with open(fqlist, 'r') as in_handle:
        for line in in_handle:
            fq_1, fq_2 = line.strip().split("\t")
            r1.append(os.path.abspath(fq_1))
            r2.append(os.path.abspath(fq_2))
    pe1 = ",".join(r1)
    pe2 = ",".join(r2)
    coasm_shell = "%s -1 %s -2 %s -t %d --out-dir %s --out-prefix %s\n" % (shutil.which("megahit"), pe1, pe2, thread, outdir, prefix)
    print(coasm_shell)

    with open("./megahit_coasm.sh", 'w') as sh_h:
        sh_h.write(coasm_shell)
    with open("./megahit_coasm_submit.sh", 'w') as sge_h:
        sge_h.write("qsub -cwd -q %s -P %s -l %s megahit_coasm.sh\n" % (queue, project, resource))

    '''
    jobname = "megahit_coasm" + "_" + datetime.now().strftime("%Y%m%d%H%M%S")
    logdir = jobname + "_qsub"
    if os.path.exists(logdir):
        os.remove(logdir)
    os.makedirs(logdir)

    jobfile = os.path.join(logdir, jobname + "_1.sh")
    with open(jobfile, 'w') as out_handle:
        out_handle.write(coasm_shell)
    
    submit_job(jobname, 1, queue, project, resource, logdir)
    '''


def main():
    parser = argparse.ArgumentParser(description='using megahit to do assembly or coassembly')
    parser.add_argument('-asm', action='store_true', help='do assembly', default=False)
    parser.add_argument('-coasm', action='store_true', help='do coassembly', default=False)
    parser.add_argument('-fqlist', type=str, help='clean pair-ended reads, each line format: reads_1.fq.gz reads_2.fq.gz')
    parser.add_argument('-thread', type=int, help="number of CPU threads, at least 2 if GPU enabled. [# of logical processors]", default=8)
    parser.add_argument('-outdir', type=str, help='output directory', default="coasm_results")
    parser.add_argument('-prefix', type=str, help='coassembly prefix', default="megahit_coasm.out")
    parser.add_argument('-queue', type=str, help='submit queue', default='st.q')
    parser.add_argument('-project', type=str, help='project id', default='F16ZQSB1SY2779')
    parser.add_argument('-resource',type=str, help='resourse requirment', default='vf=30G,p=8')

    args = parser.parse_args()

    assert int(args.resource.split("=")[2]) == args.thread, "please let p number equal thread number"

    #if args.asm:
    #    assembly(args.fqlist)
    
    if args.coasm:
        coassembly(args.fqlist, args.thread, args.outdir, args.prefix, args.queue, args.project, args.resource)

if __name__ == '__main__':
    main()