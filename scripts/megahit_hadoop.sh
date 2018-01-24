#!/bin/bash
# Author: liuxing 
# Email: liuxing2@genomics.cn

if [[ $# -ne 8 ]];then
        echo
        echo "usage: $0 -l FastaqFileList -o OutputDirPath -d HdfsOutputPath -n NumberOfTasks
        -l fastaq file list, please make a list including all fastq file, path one sample per line
           and seperate the read1 read2 and singleRead with space or table e.g: read1.fq read2.fq singleRead.fq
        -o output directory path, the directory that you would write the run script and assembly result
        -d HDFS output path, e.g: /user/liuxing2/megahitout
        -n the number of tasks, equal to the number of the samples "
        echo
else
        while [[ -n "$1" ]]
        do
                case "$1" in
                        -l) fqfilelist="$2"
                            shift ;;
                        -o) outpath="$2"
                            shift ;;
                        -d) dfsoutpath="$2"
                            shift ;;
                        -n) maps="$2"
                            shift ;;
                esac
                shift
        done
        if [[ ! -d $outpath ]];then
                mkdir $outpath
        fi
     
        echo "while read LINE
        do
                if [[ -n \$LINE ]];then
                        echo \$LINE;
                        read1=\`echo \$LINE| awk '{print \$2}'\`
                        read2=\`echo \$LINE| awk '{print \$3}'\`
			reads=\`echo \$LINE| awk '{print \$4}'\`                     
                        base=\`basename \$read1\`
                        prefix=\${base%%.*}
			outputfilename=\${prefix}.megahit_asm
                        /hwfssz1/ST_META/CD/zhujie/program/bioenv/bin/megahit -1 \$read1 -2 \$read2 -r \$reads -o ${outpath}/\$outputfilename --out-prefix \$prefix
                fi
        done" >${outpath}/megahit.sh

        echo "/hwfssz1/BIGDATA_COMPUTING/hadoop/job_submit/10.53.20.169/CDH/bin/hadoop fs -rm -r -skipTrash $dfsoutpath
/hwfssz1/BIGDATA_COMPUTING/hadoop/job_submit/10.53.20.169/CDH/bin/hadoop jar /hwfssz1/BIGDATA_COMPUTING/hadoop/job_submit/10.53.20.169/CDH/jars/hadoop-streaming-2.6.0-cdh5.11.1.jar -D mapreduce.job.name=\"megahit\" -D mapreduce.job.maps=$maps -D mapreduce.job.reduces=0  -D mapreduce.map.memory.mb=25600 -inputformat org.apache.hadoop.mapred.lib.NLineInputFormat -input file:$fqfilelist -output $dfsoutpath -mapper \"sh megahit.sh\" -file ${outpath}/megahit.sh

/hwfssz1/BIGDATA_COMPUTING/hadoop/job_submit/10.53.20.169/CDH/bin/hadoop fs -rm -r -skipTrash $dfsoutpath" >${outpath}/megahit_hadoopsubmit.sh
fi

