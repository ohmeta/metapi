#!/usr/bin/env python
'''                            _     _                                                       
                              | |   | |                                                      
   __ _ ___ ___  ___ _ __ ___ | |__ | | ___ _ __  __      ___ __ __ _ _ __  _ __   ___ _ __  
  / _` / __/ __|/ _ \ '_ ` _ \| '_ \| |/ _ \ '__| \ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__| 
 | (_| \__ \__ \  __/ | | | | | |_) | |  __/ |     \ V  V /| | | (_| | |_) | |_) |  __/ |    
  \__,_|___/___/\___|_| |_| |_|_.__/|_|\___|_|      \_/\_/ |_|  \__,_| .__/| .__/ \___|_|    
                                                                     | |   | |               
                                                                     |_|   |_|               
'''
import argparse
import os
import shutil
import subprocess
import sys
import textwrap
import pandas as pd
from argparse import RawTextHelpFormatter

__author__ = 'Jie Zhu, Xing Liu, Zhi Feng Wang'
__version__ = '2.0.0'
__date__ = 'March 30, 2019'

def parse_samples(samples_tsv):
    return pd.read_csv(samples_tsv, sep='\s+').set_index("id", drop=False)


def get_fqpath(sample_df, sample_id, col):
    return sample_df.loc[sample_id, [col]].dropna()[0]


METASPADES_TEMPLATE = '''metaspades.py \
-1 {reads1} \
-2 {reads2} \
-k {kmers} \
{only_assembler} \
--threads {threads} \
-o {out_dir} 2> {log}
pigz -p {threads} {out_dir}/scaffolds.fasta
mv {out_dir}/scaffolds.fasta.gz {scaftigs}
rm -rf {kmer_dirs}
rm -rf {out_dir}/corrected    
if {only_save_scaftigs}; then
    find {out_dir} -type f ! -wholename "{scaftigs}" -delete
else
    find {out_dir} -type f ! -wholename "{scaftigs}" ! -wholename "{tar_results}" | xargs -I % sh -c 'tar -rf {tar_results} %; rm -rf %'
    pigz -p {threads} {tar_results}
fi       
rm -rf {out_dir}/tmp
rm -rf {out_dir}/misc'''


class metaspadeser:
    def __init__(self, sample_id, reads1, reads2, only_assembler, kmers, threads, out_dir, only_save_scaftigs, log):
        self.reads1 = reads1
        self.reads2 = reads2
        self.kmers = ",".join(kmers)
        self.threads = threads
        self.out_dir = out_dir
        self.only_save_scaftigs = only_save_scaftigs
        if only_save_scaftigs:
            self.only_save_scaftigs = "true"
        else:
            self.only_save_scaftigs = "false"
        self.tar_results = os.path.join(out_dir, sample_id + ".metaspades.tar")
        self.log = log
        self.scaftigs = os.path.join(out_dir, sample_id + ".scaftigs.fa.gz")
        self.kmer_dirs = " ".join([os.path.join(out_dir, "K" + str(kmer)) for kmer in kmers])
        if only_assembler:
            self.only_assembler = "--only-assembler"
        else:
            self.only_assembler = ""


def parse_line(line, addse):
    '''parse fastq file one line'''
    if (len(line.split()) == 3) or ((len(line.split()) == 2) and (not addse)):
        return line.split()
    else:
        sys.exit(
            '''fastq path list is wrong, need two or three columns each line\n
                if addse, need three columns each line\n
                fastq_1_path\tfastq_2_path\tfastq_single_path''')


def sge_init(assembler, out_dir, job_dir, job_line, project, queue, resource):
    '''return a sge configuration dict'''
    sge = {}
    sge["asub"] = shutil.which("asub.py")
    sge["project"] = project
    sge["queue"] = queue
    sge["resource"] = resource
    sge["jobname"] = "assembly"
    sge["jobline"] = job_line
    sge["assembler"] = shutil.which(assembler)
    sge["jobdir"] = job_dir
    sge["out_dir"] = out_dir
    sge["assembly_script"] = os.path.join(sge["jobdir"],
                                          assembler + "_assembly.sh")
    sge["submit_script"] = os.path.join(sge["jobdir"],
                                        assembler + "_submit.sh")
    if sge["asub"] is None:
        sge["asub"] = os.path.join(os.path.dirname(os.path.abspath(__file__)), "asub.py")
    submit_cmd = "python %s -project %s -queue %s -jobfile %s -jobname %s -jobline %d -resource %s\n" % (
        sge["asub"], sge["project"], sge["queue"],
        os.path.basename(sge["assembly_script"]), sge["jobname"],
        sge["jobline"], sge["resource"])
    with open(sge["submit_script"], 'w') as submit_h:
        submit_h.write(submit_cmd)
    return sge


def gen_metaspades_sge(sge, samples, only_assembler, klist, only_save_scaftigs, threads, logdir):
    samples_df = parse_samples(samples)
    with open(sge["assembly_script"], 'w') as assembly_handle:
        for sample_id in samples_df.index:
            r1 = get_fqpath(samples_df, sample_id, "fq1")
            r2 = get_fqpath(samples_df, sample_id, "fq2")
            out_dir_asm = os.path.join(sge["out_dir"],
                                       sample_id + ".metaspades_out")
            log = os.path.join(logdir, sample_id + ".metaspades.log")
            metaspades_cmd = METASPADES_TEMPLATE.format_map(
                vars(metaspadeser(sample_id, r1, r2, only_assembler, klist, threads, out_dir_asm, only_save_scaftigs, log)))
            assembly_handle.write(metaspades_cmd + "\n")


def gen_megahit_sge(sge, addse, fq_list, nomercy, mincount, klist, kmin, kmax,
                    kstep, threads):
    '''generate megahit script for sge'''
    with open(sge["assembly_script"], 'w') as assembly_handle:
        with open(fq_list, 'r') as fqlist_handle:
            for line in fqlist_handle:
                fq_path = parse_line(line.strip(), addse)
                sample_name = os.path.basename(fq_path[1]).split('.')[0]
                out_dir_asm = os.path.join(sge["out_dir"],
                                           sample_name + ".megahit_out")
                assembly_shell = sge["assembler"] + " --min-count " + str(
                    mincount)
                if not klist:
                    assembly_shell += " --k-min " + str(kmin) + \
                        " --k-max " + str(kmax) + \
                        " --k-step " + str(kstep)
                else:
                    assembly_shell += " --k-list " + klist
                if nomercy:
                    assembly_shell += " --no-mercy "
                assembly_shell += " -1 " + fq_path[0] + \
                                  " -2 " + fq_path[1]
                if addse:
                    assembly_shell += " -r " + fq_path[2]
                assembly_shell += " -t " + str(threads) + \
                                  " --out-dir " + out_dir_asm + \
                                  " --out-prefix " + sample_name + '\n'
                assembly_handle.write(assembly_shell)


def gen_idba_ud_sge(sge, addse, precor, fq_list, kmin, kmax, kstep, mincontig,
                    threads):
    '''generate idba_ud script for sge'''
    fq2fa = shutil.which("fq2fa")
    with open(sge["assembly_script"], 'w') as assembly_handle:
        with open(fq_list, 'r') as fqlist_handle:
            for line in fqlist_handle:
                fq_path = parse_line(line.strip(), addse)
                # print(fq_path)
                sample_name = os.path.basename(fq_path[1]).split('.')[0]
                fq_1 = sample_name + ".1.fq"
                fq_2 = sample_name + ".2.fq"
                fa_pe = sample_name + ".pe.fa"
                fq_se = sample_name + ".se.fq"
                fa_se = sample_name + ".se.fa"
                out_dir_asm = os.path.join(sge["out_dir"],
                                           sample_name + ".idba_ud_out")
                assembly_shell = "mkdir " + out_dir_asm + "\n" + \
                                 "cd " + out_dir_asm + "\n" + \
                                 "zcat " + fq_path[0] + " > " + fq_1 + "\n" + \
                                 "zcat " + fq_path[1] + " > " + fq_2 + "\n" + \
                                 fq2fa + " --merge " + fq_1 + " " + fq_2 + " " + fa_pe + "\n"
                if addse:
                    assembly_shell += "zcat " + fq_path[2] + " > " + fq_se + "\n" + \
                                      fq2fa + " " + fq_se + " " + fa_se + "\n" + \
                                      "rm -f *.fq\n" + \
                                      sge["assembler"] + " -r " + fa_pe + " -l " + fa_se + \
                                      " --mink " + str(kmin) + " --maxk " + str(kmax) + " --step " + str(kstep) + \
                                      " --min_contig " + str(mincontig) + \
                                      " -o " + "./" + \
                                      " --num_threads " + str(threads)
                else:
                    assembly_shell += "rm -f *.fq\n" + \
                                      sge["assembler"] + " -r " + fa_pe + \
                                      " --mink " + str(kmin) + " --maxk " + str(kmax) + " --step " + str(kstep) + \
                                      " --min_contig " + str(mincontig) + \
                                      " -o " + "./" + \
                                      " --num_threads " + str(threads)
                if precor:
                    assembly_shell += " --pre_correction "
                if addse:
                    assembly_shell += "\nrm -r kmer contig-* align-* graph-* local-contig-* " + \
                                      fa_pe + " " + fa_se + "\n"
                else:
                    assembly_shell += "\nrm -r kmer contig-* align-* graph-* local-contig-* " + \
                                      fa_pe + "\n"
                final_asm_fa = sample_name + ".scaffold.fa.gz"
                assembly_shell += "pigz scaffold.fa\nmv scaffold.fa.gz %s\n" % (
                    final_asm_fa)

                assembly_handle.write(assembly_shell)


def gen_megahit_hadoop(assembler, addse, fq_list, out_dir, mincount, kmin,
                       kmax, kstep, threads, job_dir, memory):
    '''generate megahit script for hadoop'''
    # fq_list = os.path.abspath(fq_list)
    # job_dir = os.path.abspath(job_dir)
    # out_dir = os.path.abspath(out_dir)
    megahit = shutil.which('megahit')
    handle = open(fq_list, 'r')
    sample_num = len(handle.readlines())
    handle.close()

    submit_script = os.path.join(job_dir, assembler + "_hadoop_submit.sh")
    assembly_script = os.path.join(job_dir,
                                   assembler + "_assembly_template.sh")

    hdfs_outdir = os.path.join(job_dir, "hdfs_temp")
    if os.path.exists(hdfs_outdir):
        os.removedirs(hdfs_outdir)

    hadoop = {}
    hadoop["bin"] = shutil.which("hadoop")
    hadoop["jar"] = shutil.which("hadoop-streaming-2.6.0-cdh5.11.1.jar")
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
        if addse:
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
        %s --min-count %s --k-min %s --k-max %s --k-step %s -1 $read1 -2 $read2 -r $reads -t %s --out-dir %s/$outputfilename --out-prefix $prefix
    fi
done\n''' % (megahit, mincount, kmin, kmax, kstep, threads, out_dir)
        else:
            assembly_shell = '''while read LINE
do
    if [[ -n $LINE ]];then
        echo $LINE;
        read1=`echo $LINE| awk '{print $2}'`
        read2=`echo $LINE| awk '{print $3}'`
        base=`basename $read1`
        prefix=${base%%%%.*}
        outputfilename=${prefix}.megahit_out
        %s --min-count %s --k-min %s --k-max %s --k-step %s -1 $read1 -2 $read2 -t %s --out-dir %s/$outputfilename --out-prefix $prefix
    fi
done\n''' % (megahit, mincount, kmin, kmax, kstep, threads, out_dir)

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


def main():
    '''main function'''
    banner = '''
            ================================================================================================================
            Yet another genome assembler wrapper for PE + SE reads assembly using megahit or idba_ud based on SGE or Hadoop

            GPL License
            ================================================================================================================
'''
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description=textwrap.dedent(banner))
    parser.add_argument(
        '--assembler',
        type=str,
        default='megahit',
        help='assembler supported, megahit or idba_ud, default: megahit')
    parser.add_argument(
        '--platform',
        type=str,
        default='sge',
        help='platform, sge or hadoop, default: sge')
    parser.add_argument(
        '--fqlist',
        type=str,
        help='rmhost fastq file list, two or three columns each line')
    parser.add_argument(
        '--addse',
        action='store_true',
        help=
        'add single fastq, if --addse: yes(fqlist need three columns each line), else: no'
    )
    parser.add_argument(
        '--kmin',
        type=int,
        default=21,
        help='kmer min length, must be odd number, default: 21')
    parser.add_argument(
        '--kmax',
        type=int,
        default=99,
        help='kmer max length, must be odd number, default: 99')
    parser.add_argument(
        '--kstep',
        type=int,
        default=10,
        help='kmer step, must be even number, default: 10')
    parser.add_argument(
        '--klist',
        type=str,
        nargs='*',
        help='all must be odd, in the range 15-255, increment <= 28, [21,29,39,59,79,99,119,141]'
    )

    metaspades_group = parser.add_argument_group("metaspades", "args for metaspades")
    metaspades_group.add_argument(
        '--samples',
        type=str,
        help='samples tsv, id, fq1, fq2'
    )
    metaspades_group.add_argument(
        '--only_assembler',
        default=False,
        action='store_true',
        help='runs only assembling(without read error correction)'
    )
    metaspades_group.add_argument(
        '--only_save_scaftigs',
        default=False,
        action='store_true',
        help='only save scaftigs'
    )
    metaspades_group.add_argument(
        '--asmlogdir',
        type=str,
        help='metaspades assembly log dir'
    )

    megahit_group = parser.add_argument_group("megahit", "args for megahit")
    megahit_group.add_argument(
        '--mincount',
        type=int,
        default=2,
        help='minimum multiplicity for filtering (k_min+1)-mers'
    )
    megahit_group.add_argument(
        '--nomercy',
        action='store_true',
        help=
        'for generic dataset >= 30x, MEGAHIT may generate better results with --no-mercy option'
    )

    idba_ud_group = parser.add_argument_group("idba_ud", "args for idba_ud")
    idba_ud_group.add_argument("--pre_correction", action='store_true')
    idba_ud_group.add_argument(
        "--min_contig",
        type=int,
        help='minimum size of contig, default: 200',
        default=200)

    sge_group = parser.add_argument_group("sge", "args for sge")
    sge_group.add_argument('--project_name', type=str, help='set project name')
    sge_group.add_argument(
        '--queue',
        type=str,
        default='st.q',
        help='bind job to queue(s), default: st.q')
    sge_group.add_argument(
        '--resource_list',
        type=str,
        default='vf=10G,p=8',
        help='request the given resources, default: vf=10G,p=8')
    sge_group.add_argument(
        '--threads',
        type=int,
        default=8,
        help='number of CPU threads, for megahit at least 2 if GPU enabled')
    sge_group.add_argument('--outdir', type=str, help='assembly result outdir')
    sge_group.add_argument(
        '--jobdir',
        type=str,
        default='job',
        help='a dir contain assembly script and submit script, default: job')

    hadoop_group = parser.add_argument_group("hadoop", "args for hadoop")
    hadoop_group.add_argument(
        '--memory',
        type=int,
        default=25000,
        help='set memory(MB) usage for hadoop job, default: 25000')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    os.makedirs(args.jobdir, exist_ok=True)

    if args.platform == "sge":
        if int(args.resource_list.split("=")[-1]) != int(args.threads):
            sys.exit("Please let sge resource_list p = threads")

    if args.platform == "sge" and args.assembler == "metaspades":
        sge = sge_init(args.assembler, args.outdir, args.jobdir, 13,
                       args.project_name, args.queue, args.resource_list)
        gen_metaspades_sge(sge, args.samples, args.only_assembler,
                           args.klist, args.only_save_scaftigs, args.threads, args.asmlogdir)

    elif args.platform == "sge" and args.assembler == "megahit":
        sge = sge_init(args.assembler, args.outdir, args.jobdir, 1,
                       args.project_name, args.queue, args.resource_list)
        gen_megahit_sge(sge, args.addse, args.fqlist, args.nomercy,
                        args.mincount, args.klist, args.kmin, args.kmax,
                        args.kstep, args.threads)

    elif args.platform == "sge" and args.assembler == "idba_ud":
        sge = {}
        if args.addse:
            sge = sge_init(args.assembler, args.outdir, args.jobdir, 11,
                           args.project_name, args.queue, args.resource_list)
        else:
            sge = sge_init(args.assembler, args.outdir, args.jobdir, 10,
                           args.project_name, args.queue, args.resource_list)
        gen_idba_ud_sge(sge, args.addse, args.pre_correction, args.fqlist,
                        args.kmin, args.kmax, args.kstep, args.min_contig,
                        args.threads)

    elif args.platform == "hadoop" and args.assembler == "megahit":
        gen_megahit_hadoop(args.assembler, args.addse, args.fqlist,
                           args.outdir, args.mincount, args.kmin, args.kmax,
                           args.kstep, args.threads, args.jobdir, args.memory)
    else:
        print(
            args.platform + " and " + args.assembler +
            " have not been implemented\nThanks using assembler_wrapper.py\n")


if __name__ == "__main__":
    main()
'''                            _     _                                                       
                              | |   | |                                                      
   __ _ ___ ___  ___ _ __ ___ | |__ | | ___ _ __  __      ___ __ __ _ _ __  _ __   ___ _ __  
  / _` / __/ __|/ _ \ '_ ` _ \| '_ \| |/ _ \ '__| \ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__| 
 | (_| \__ \__ \  __/ | | | | | |_) | |  __/ |     \ V  V /| | | (_| | |_) | |_) |  __/ |    
  \__,_|___/___/\___|_| |_| |_|_.__/|_|\___|_|      \_/\_/ |_|  \__,_| .__/| .__/ \___|_|    
                                                                     | |   | |               
                                                                     |_|   |_|               
 Run assembler wrapper for all samples http://patorjk.com/software/taag
'''
