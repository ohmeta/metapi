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

__author__ = 'Jie Zhu'
__version__ = '0.1.0'
__date__ = 'March 24, 2018'


def parse_job(job_name, job_file):
    # how to remove job_file parameter
    with fileinput.input(files=job_file if not job_file is None else ('-', )) as in_h:
        if fileinput.filename() is None:
            job_f = os.path.join(os.path.curdir, job_name.rstrip(".sh") + ".sh")
            with open(job_f, 'w') as job_h:
                for one_line in in_h:
                    #if one_line.strip() and (not one_line.startswith("#")):
                    job_h.write(one_line)
            total_line = fileinput.lineno()
            return job_f, total_line
        else:
            total_line = len(open(fileinput.filename(), 'r').readlines())
            return fileinput.filename(), total_line


# method_1
# -t 1-$total_job_line:$a_job_line
# begin=$SGE_TASK_ID
# end=$begin+$a_job_line-1

# method_2 (choose this method)
# -t 1-($total_job_line/$a_job_line):1
# end=$SGE_TASK_ID*a_job_line
# begin=$end-$a_job_line+1


def submit_job(job_f, job_name, total_job_line, a_job_line, queue, prj_id, resource, logdir):
    if total_job_line == 0:
        print("nothing to do!")
        sys.exit(-1)
    elif total_job_line % a_job_line == 0:
        array_size = total_job_line / a_job_line
        array_range = "1-" + str(int(array_size)) + ":1"
    elif total_job_line < a_job_line:
        array_range = "1-1:1"
        a_job_line = total_job_line
    else:
        print("wrong [total job line: %s] and [a job line: %s]" % (total_job_line, a_job_line))
        sys.exit(-1)

    submit_f = os.path.join(os.path.curdir, job_name.rstrip(".sh") + "_submit.sh")
    with open(submit_f, 'w') as submit_h:
        submit_h.write('''#!/bin/bash\n\
#$ -S /bin/bash
#$ -N %s
#$ -cwd
#$ -l %s
#$ -q %s
#$ -P %s
#$ -t %s
jobfile=%s
((end=$SGE_TASK_ID*%d))
((begin=$end-%d+1))
cmd=$(sed -n -e \"$begin,$end p\" $jobfile)
$cmd\n''' % (job_name, resource, queue, prj_id, array_range, job_f, a_job_line, a_job_line))

    os.chmod(submit_f, stat.S_IRWXU)
    submit_cmd = shutil.which("qsub") + \
                 " -e " + logdir + job_name + "_\\$TASK_ID.e" + \
                 " -o " + logdir + job_name + "_\\$TASK_ID.o " + submit_f
    print(submit_cmd)
    subprocess.call(submit_cmd, shell=True)

def main():
    '''it is a very simple script to submit array job, but you need supply real run command'''
    parser = argparse.ArgumentParser(description='make submit array job easy')
    # how to remove -jobfile parameter
    parser.add_argument('-jobfile', nargs='*', help='job file to read, if empty, stdin is used')
    parser.add_argument('-jobname', type=str, help='job name', default='job')
    parser.add_argument('-jobline', type=int, help='set the number of lines to form a job', default=1)
    parser.add_argument('-queue', type=str, help='submit queue', default='st.q')
    parser.add_argument('-project', type=str, help='project id', default='F16ZQSB1SY2779')
    parser.add_argument('-resource', type=str, help='resourse requirment', default='vf=128M,p=1')
    parser.add_argument('-logdir', type=str, help='array job log directory')
    args = parser.parse_args()

    assert re.match(r'vf=\d+\w,p=\d+', args.resource), "please specific memory usage and number processor"
    assert not re.match(r'^\d+', args.jobname), "array job name cannot start with a digit"
    assert args.jobline >= 1, "a job line can't to be zero"
    
    args.jobname += "_" + datetime.now().strftime("%Y%m%d%H%M%S")
    
    if not args.logdir:
        args.logdir = args.jobname + "_qsub"
    args.logdir = args.logdir.rstrip("/") + "/"
    
    if os.path.exists(args.logdir):
        os.remove(args.logdir)
    os.makedirs(args.logdir)

    job_f, total_job_line = parse_job(args.jobname, args.jobfile)
    submit_job(job_f, args.jobname, total_job_line, args.jobline, args.queue, args.project, args.resource, args.logdir)


if __name__ == '__main__':
    main()
