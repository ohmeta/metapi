#!/usr/bin/env python
'''
  __  __ ______ _____          _    _ _____ _______  __          _______            _____  _____  ______ _____
 |  \/  |  ____/ ____|   /\   | |  | |_   _|__   __| \ \        / /  __ \     /\   |  __ \|  __ \|  ____|  __ \
 | \  / | |__ | |  __   /  \  | |__| | | |    | |     \ \  /\  / /| |__) |   /  \  | |__) | |__) | |__  | |__) |
 | |\/| |  __|| | |_ | / /\ \ |  __  | | |    | |      \ \/  \/ / |  _  /   / /\ \ |  ___/|  ___/|  __| |  _  /
 | |  | | |___| |__| |/ ____ \| |  | |_| |_   | |       \  /\  /  | | \ \  / ____ \| |    | |    | |____| | \ \
 |_|  |_|______\_____/_/    \_\_|  |_|_____|  |_|        \/  \/   |_|  \_\/_/    \_\_|    |_|    |______|_|  \_\

https://github.com/voutcn/megahit/wiki/Assembly-Tips

Choosing k:
    MEGAHIT uses multiple k-mer strategy.
    Minimum k, maximum k and the step for iteration can be set by options
    --k-min, --k-max and --k-step respectively.
    k must be odd numbers while the step must be an even number.
    - for ultra complex metagenomics data such as soil, a larger kmin, say 27,
      is recommended to reduce the complexity of the de Brujin graph.
      Quality trimming is also recommended.
    - for high-depth generic data, large --k-min (25 to 31) is recommended.
    - smaller --k-step, say 10, is more friendly to low-coverage datasets.

Filtering (kmin+1)-mer:
    (kmin+1)-mer with multiplicity lower than d (default 2, specified by --min-count option)
    will be discarded. You should be cautious to set d less than 2,
    which will lead to a much larger and noisy graph.
    We recommend using the default value 2 for metagenomics assembly.
    If you want to use MEGAHIT to do generic assemblies, please change this value
    according to the sequencing depth(recommend --min-count 3 for >40x).

Mercy k-mer:
    This is specially designed for metagenomics assembly to recover low coverage sequence.
    For generic dataset >= 30x, MEGAHIT may generate better results with --no-mercy option.

k-min 1pass mode:
    This mode can be activated by option --kmin-1pass.
    It is more memory efficient for ultra low-depth datasets, such as soil metagenomics data.

k must be odd numbers while the step must be an even number
    default(bp):
        k_min = 21
        k_max = 141
        k_step = 10
        k_list = [21,29,39,59,79,99,119,141]
        min_count = 2
        min_contig_len = 200
    meta-sensitive(bp):
        min_count = 1
        k_list = [21,29,39,49,59,69,79,89,99,109,119,129,141]
        set_list_by_min_max_step = False
    meta-large(bp):
        k_min = 27
        k_max = 127
        k_step = 10
        min_count = 1
        set_list_by_min_max_step = True
'''

import argparse
from argparse import RawTextHelpFormatter
import textwrap
import os
import stat
import shutil

__author__ = 'Jie Zhu (zhujie@genomics.cn), Xing Liu (lingxing2@genomics.cn)'
__version__ = '1.0.0'
__date__ = 'January 16, 2018'


def gen_megahit_sge(fq_list, out_dir, kmin, kmax, kstep,
                    job_dir, project_name, queue, resource_list):
    '''generate megahit script for sge'''
    megahit = shutil.which('megahit')
    out_dir = os.path.abspath(out_dir)
    job_dir = os.path.abspath(job_dir)
    submit_script = os.path.join(job_dir, "megahit_sge_submit.sh")
    assembly_script = os.path.join(job_dir, "megahit_assembly.sh")

    sge = {}
    sge["qsub_all_path"] = "/hwfssz1/ST_META/CD/zhujie/bin/qsub_all.pl"
    sge["project_name"] = project_name
    sge["queue_name"] = queue
    sge["resource_list"] = resource_list
    sge["log_dir"] = os.path.join(job_dir, "qsub_log")
    sge["job_name"] = "megahit"
    sge["job_number"] = str(80)
    sge["script_path"] = assembly_script

    with open(assembly_script, 'w') as assembly_handle:
        with open(fq_list, 'r') as fqlist_handle:
            for line in fqlist_handle:
                (reads_a, reads_b, reads_s) = line.split()
                reads_a = os.path.abspath(reads_a)
                reads_b = os.path.abspath(reads_b)
                reads_s = os.path.abspath(reads_s)
                sample_name = os.path.basename(reads_a).split('.')[0]
                out_dir_asm = os.path.join(
                    out_dir, sample_name + ".megahit_out")
                assembly_shell = megahit + \
                    " --k-min " + str(kmin) + \
                    " --k-max " + str(kmax) + \
                    " --k-step " + str(kstep) + \
                    " -1 " + reads_a + \
                    " -2 " + reads_b + \
                    " -r " + reads_s + \
                    " --outdir " + out_dir_asm + \
                    " --out-prefix " + sample_name + '\n'
                assembly_handle.write(assembly_shell)

    with open(submit_script, 'w') as submit_handle:
        submit_shell = "perl " + sge["qsub_all_path"] + \
                       " -q " + sge["queue_name"] + \
                       " -P " + sge["project_name"] + \
                       " -l " + sge["resource_list"] + \
                       " -N " + sge["job_name"] + \
                       " -m " + sge["job_number"] + \
                       " -d " + sge["log_dir"] + \
                       " " + sge["script_path"] + '\n'
        submit_handle.write(submit_shell)

    os.chmod(assembly_script, stat.S_IRWXU)
    os.chmod(submit_script, stat.S_IRWXU)


def gen_megahit_hadoop(fq_list, out_dir, kmin, kmax, kstep, job_dir, memory):
    '''generate megahit script for hadoop'''
    fq_list = os.path.abspath(fq_list)
    job_dir = os.path.abspath(job_dir)
    out_dir = os.path.abspath(out_dir)
    megahit = shutil.which('megahit')
    handle = open(fq_list, 'r')
    sample_num = len(handle.readlines())
    handle.close()

    submit_script = os.path.join(job_dir, "megahit_hadoop_submit.sh")
    assembly_script = os.path.join(job_dir, "megahit_assembly_template.sh")

    hdfs_outdir = os.path.join(job_dir, "hdfs_temp")
    if os.path.exists(hdfs_outdir):
        os.removedirs(hdfs_outdir)

    hadoop = {}
    hadoop["bin"] = "/hwfssz1/BIGDATA_COMPUTING/hadoop/job_submit/10.53.20.169/CDH/bin/hadoop"
    hadoop["jar"] = "/hwfssz1/BIGDATA_COMPUTING/hadoop/job_submit/10.53.20.169/CDH/jars/hadoop-streaming-2.6.0-cdh5.11.1.jar"
    hadoop["job_name"] = "\"megahit_assembly\""
    hadoop["maps"] = str(sample_num)
    hadoop["reduces"] = str(0)
    hadoop["memory"] = str(memory)
    hadoop["input_format"] = "org.apache.hadoop.mapred.lib.NLineInputFormat"
    hadoop["input"] = "file:" + fq_list
    hadoop["output"] = "file:" + hdfs_outdir
    hadoop["mapper"] = "\"bash megahit_assembly.sh\""
    hadoop["file"] = assembly_script

    with open(assembly_script, 'w') as assembly_handle:
        assembly_shell = '''while read LINE
do
    if [[ -n $LINE ]];then
        echo $LINE;
        read1=`echo $LINE| awk '{print $2}'`
        read2=`echo $LINE| awk '{print $3}'`
        reads=`echo $LINE| awk '{print $4}'`
        base=`basename $read1`
        prefix=${base%%%%.*}
        outputfilename=${prefix}.megahit_out
        %s --k-min %s --k-max %s --k-step %s -1 $read1 -2 $read2 -r $reads --out-dir %s/$outputfilename --out-prefix $prefix
    fi
done\n''' % (megahit, kmin, kmax, kstep, out_dir)
        assembly_handle.write(assembly_shell)

    with open(submit_script, 'w') as submit_handle:
        submit_shell = hadoop["bin"] + \
            " jar " + hadoop["jar"] + \
            " -D mapreduce.job.name=" + hadoop["job_name"] + \
            " -D mapreduce.job.maps=" + hadoop["maps"] + \
            " -D mapreduce.job.reduces=" + hadoop["reduces"] + \
            " -D mapreduce.map.memory.mb=" + hadoop["memory"] + \
            " -inputformat " + hadoop["input_format"] + \
            " -input " + hadoop["input"] + \
            " -output " + hadoop["output"] + \
            " -mapper " + hadoop["mapper"] + \
            " -file " + hadoop["file"] + '\n'
        submit_handle.write(submit_shell)
    os.chmod(assembly_script, stat.S_IRWXU)
    os.chmod(submit_script, stat.S_IRWXU)


def main():
    '''main function'''
    banner = '''
            ===========================================================================================
            Yet another megahit wrapper for PE + SE reads assembly using megahit based on SGE or Hadoop

            GPL License

            zhujie@genomics.cn liuxing2@genomics.cn
            ===========================================================================================
'''
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description=textwrap.dedent(banner))

    parser.add_argument('--platform', type=str, default='sge',
                        help='platform, sge or hadoop, default: sge')
    parser.add_argument('--fqlist', type=str,
                        help='rmhost fastq file list')
    parser.add_argument('--kmin', type=int, default=21,
                        help='kmer min length, must be odd number, default: 21')
    parser.add_argument('--kmax', type=int, default=99,
                        help='kmer max length, must be odd number, default: 99')
    parser.add_argument('--kstep', type=int, default=10,
                        help='kmer step, must be even number, default: 10')
    parser.add_argument('--outdir', type=str, default='megahit_out',
                        help='assembly result outdir, default: megahit_out')
    parser.add_argument('--jobdir', type=str, default='job',
                        help='a dir contain assembly script and submit script, default: job')

    sge_group = parser.add_argument_group("sge", "args for sge")
    sge_group.add_argument('--project_name', type=str, default='F16ZQSB1SY2779',
                           help='set job\'s project, default: F16ZQSB1SY2779')
    sge_group.add_argument('--queue', type=str, default='st.q',
                           help='bind job to queue(s), default: st.q')
    sge_group.add_argument('--resource_list', type=str, default='vf=10G,p=8',
                           help='request the given resources, default: vf=10G,p=8')

    hadoop_group = parser.add_argument_group("hadoop", "args for hadoop")
    hadoop_group.add_argument('--memory', type=int, default=25000,
                              help='set memory(MB) usage for hadoop job, default: 25000')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    os.makedirs(args.jobdir, exist_ok=True)

    if args.platform == "sge":
        gen_megahit_sge(args.fqlist, args.outdir,
                        args.kmin, args.kmax, args.kstep,
                        args.jobdir, args.project_name, args.queue, args.resource_list)

    if args.platform == "hadoop":
        gen_megahit_hadoop(args.fqlist, args.outdir,
                           args.kmin, args.kmax, args.kstep, args.jobdir, args.memory)


if __name__ == "__main__":
    main()
'''
  __  __ ______ _____          _    _ _____ _______  __          _______            _____  _____  ______ _____
 |  \/  |  ____/ ____|   /\   | |  | |_   _|__   __| \ \        / /  __ \     /\   |  __ \|  __ \|  ____|  __ \
 | \  / | |__ | |  __   /  \  | |__| | | |    | |     \ \  /\  / /| |__) |   /  \  | |__) | |__) | |__  | |__) |
 | |\/| |  __|| | |_ | / /\ \ |  __  | | |    | |      \ \/  \/ / |  _  /   / /\ \ |  ___/|  ___/|  __| |  _  /
 | |  | | |___| |__| |/ ____ \| |  | |_| |_   | |       \  /\  /  | | \ \  / ____ \| |    | |    | |____| | \ \
 |_|  |_|______\_____/_/    \_\_|  |_|_____|  |_|        \/  \/   |_|  \_\/_/    \_\_|    |_|    |______|_|  \_\

 Run megahit wrapper for all samples http://patorjk.com/software/taag
'''
