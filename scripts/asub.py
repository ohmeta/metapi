#!/usr/bin/env python
# please see https://github.com/lh3/asub
import argparse
import fileinput
import os
import re
import shutil
import stat
import subprocess
import sys
from datetime import datetime

__author__ = 'Jie Zhu, Jiahui Zhu'
__email__ = 'zhujie@genomics.cn, zhujiahui@genomics.cn'
__version__ = '0.3.1'
__date__ = 'Jun 19, 2018'


def parse_job(job_name, job_file, a_job_line, logdir):
    with fileinput.input(files=job_file if not job_file is None else ('-', )) as in_h:
        job_num = 0
        for one_line in in_h:
            job_num += 1
            job_f = os.path.join(logdir, job_name.rstrip(".sh") + "_" + str(job_num) + ".sh")
            with open(job_f, 'w') as job_h:
                job_h.write(one_line)
                while fileinput.lineno() % a_job_line != 0:
                    job_h.write(next(in_h))
                #for i in range(1, a_job_line):
                #    job_h.write(next(in_h))
        return job_num


def submit_job(job_name, total_job_num, queue, prj_id, resource, logdir):
    submit_f = os.path.join(os.path.curdir, job_name.rstrip(".sh") + "_submit.sh")
    array_range = "1-" + str(total_job_num) + ":1"
    job_script = os.path.join(logdir, job_name.rstrip(".sh") + "_$SGE_TASK_ID.sh")
    num_proc = resource.split('=')[-1]
    with open(submit_f, 'w') as submit_h:
        submit_h.write('''#!/bin/bash\n\
#$ -clear
#$ -S /bin/bash
#$ -N %s
#$ -cwd
#$ -l %s
#$ -binding linear:%s
#$ -q %s
#$ -P %s
#$ -t %s
jobscript=%s
bash $jobscript\n''' % (job_name, resource, num_proc, queue, prj_id, array_range, job_script))

    os.chmod(submit_f, stat.S_IRWXU)
    submit_cmd = shutil.which("qsub") + \
                 " -e " + os.path.join(logdir, job_name + "_\\$TASK_ID.e") + \
                 " -o " + os.path.join(logdir, job_name + "_\\$TASK_ID.o") + " " + submit_f
    print(submit_cmd)
    subprocess.call(submit_cmd, shell=True)

def main():
    '''it is a very simple script to submit array job, but you need supply real run command'''
    parser = argparse.ArgumentParser(description='make submit array job easy')
    parser.add_argument('-jobfile', nargs='*', help='job file to read, if empty, stdin is used')
    parser.add_argument('-jobname', type=str, help='job name', default='job')
    parser.add_argument('-jobline', type=int, help='set the number of lines to form a job', default=1)
    parser.add_argument('-queue', type=str, help='submit queue', default='st.q')
    parser.add_argument('-project', type=str, help='project id', default='F16ZQSB1SY2779')
    parser.add_argument('-resource',type=str, help='resourse requirment', default='vf=1G,p=1')
    parser.add_argument('-logdir', type=str, help='array job log directory')
    args = parser.parse_args()

    assert re.match(r'vf=[\d\.]+\w,p=\d+', args.resource), "please specific memory usage and number processor"
    assert not re.match(r'^\d+', args.jobname), "array job name cannot start with a digit"
    assert args.jobline >= 1, "a job line can't to be zero"

    args.jobname += "_" + datetime.now().strftime("%Y%m%d%H%M%S")

    if not args.logdir:
        args.logdir = args.jobname + "_qsub"
    args.logdir = args.logdir.rstrip("/") + "/"

    if os.path.exists(args.logdir):
        os.remove(args.logdir)
    os.makedirs(args.logdir)

    total_job_num = parse_job(args.jobname, args.jobfile, args.jobline, args.logdir)
    submit_job(args.jobname, total_job_num, args.queue, args.project, args.resource, args.logdir)


if __name__ == '__main__':
    main()
