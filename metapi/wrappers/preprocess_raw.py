#!/usr/bin/env python

import os
import filecmp
import gzip
import subprocess as sp


input_fq_list = snakemake["input"]
output_fq_list = snakemake["output"]["reads"]

threads = snakemake["threads"]
log = snakemake["log"]

output_dir = snakemake["params"]["output_dir"]
reads_format = snakeamke["params"]["reads_format"]
is_pe = snakemake["params"]["is_pe"]
is_interleaved = snakemake["params"]["is_interleaved"]

reads_num = len(input_fq_list)

if reads_format == "fastq":
    if is_pe:
        r1 = output_fq_list[0]
        r2 = output_fq_list[1]
        fq1_ = f'''{r1}.temp.gz'''
        fq2_ = f'''{r2}.temp.gz'''
        fq_1 = f'''{r1}.paired.gz'''
        fq_2 = f'''{r2}.paired.gz'''

        if not is_interleaved:
            if reads_num == 2:
                fq1 = os.path.realpath(input_fq_list[0])
                fq2 = os.path.realpath(input_fq_list[1])
                sp.run(f'''ln -s {fq1} {fq1_} 2>{log}''', shell=True)
                sp.run(f'''ln -s {fq2} {fq2_} 2>>{log}''', shell=True)
            else:
                fq1_str = " ".join(input_fq_list[0:reads_num//2])
                fq2_str = " ".join(input_fq_list[reads_num//2:])
                sp.run(f'''cat {fq1_str} > {fq1_} 2> {log}''', shell=True)
                sp.run(f'''cat {fq2_str} > {fq2_} 2>> {log}''', shell=True)

            id1 = f'''{output_dir}/id.list.1'''
            id2 = f'''{output_dir}/id.list.2'''
            idp = f'''{output_dir}/id.list.paired'''

            sp.run(f'''seqkit seq -ni {fq1_} | sed 's#/1$##g' > {id1} 2>> {log}''', shell=True) 
            sp.run(f'''seqkit seq -ni {fq2_} | sed 's#/2$##g' > {id2} 2>> {log}''', shell=True) 

            if filecmp.cmp(id1, id2): 
                sp.run(f'''mv {fq1_} {r1} 2>> {log}''', shell=True)
                sp.run(f'''mv {fq2_} {r2} 2>> {log}''', shell=True)

                sp.run(f'''rm -rf {id1} 2>> {log}''', shell=True)
                sp.run(f'''rm -rf {id2} 2>> {log}''', shell=True)
            else:
                sp.run(
                    f'''
                    sort -T {output_dir} {id1} {id2} | \
                    uniq -c | \
                    awk '$1==2{{print $2}}' > {idp} 2>> {log}
                    ''', shell=True)

                oneline = gzip.open(fq1_, 'rt').readline().strip().split()[0]
                if "/1" in oneline:
                    sp.run(
                        f'''seqkit grep -f <(awk '{{print $0 "/1"}}' {idp} {fq1_} -o {fq_1} 2>> {log}''',
                        shell=True)
                    sp.run(
                        f'''seqkit grep -f <(awk '{{print $0 "/2"}}' {idp} {fq2_} -o {fq_2} 2>> {log}''',
                        shell=True)
                else:
                    sp.run(f'''seqkit grep -f {idp} {fq1_} -o {fq_1} 2>> {log}''', shell=True)
                    sp.run(f'''seqkit grep -f {idp} {fq2_} -o {fq_2} 2>> {log}''', shell=True)

                sp.run(f'''rm -rf {fq1_} 2>> {log}''', shell=True)
                sp.run(f'''rm -rf {fq2_} 2>> {log}''', shell=True)
                sp.run(f'''rm -rf {idp} 2>> {log}''', shell=True)
                sp.run(f'''rm -rf {id1} 2>> {log}''', shell=True)
                sp.run(f'''rm -rf {id2} 2>> {log}''', shell=True)

                ## more check
                sp.run(f'''mkdir -p {output_dir}/tmpfq''', shell=True)
                sp.run(
                    f'''
                    pigz -b 102400 -dc {fq_1} | \
                    paste - - - - | \
                    LC_ALL=C sort --parallel 4 -n -T {output_dir}/tmpfq -S 8G | \
                    tr '\t' '\n' | \
                    pigz -b 102400 > {r1} \
                    2>> {log}''', shell=True)
                sp.run(
                    f'''
                    pigz -b 102400 -dc {fq_2} | \
                    paste - - - - | \
                    LC_ALL=C sort --parallel 4 -n -T {output_dir}/tmpfq -S 8G | \
                    tr '\t' '\n' | \
                    pigz -b 102400 > {r2} \
                    2>> {log}''', shell=True)

                sp.run(f'''rm -rf {fq_1} 2>> {log}''', shell=True)
                sp.run(f'''rm -rf {fq_2} 2>> {log}''', shell=True)

        else:
            sp.run(
                f'''
                cat {input_fq_list[0]} | \
                tee >(seqtk seq -1 - | pigz -cf -p {threads} > {r1}) | \
                seqtk seq -2 - | pigz -cf -p {threads} > {r2} 2>> {log}
                ''', shell=True)
    else:
        fq_ = output_fq_list[0]
        if reads_num == 1:
            fq = os.path.realpath(input_fq_list[0])
            sp.run(f'''ln -s {fq} {fq_} 2>> {log}''', shell=True)
        else:
            fq_str = " ".join(input_fq_list)
            sp.run(f'''cat {fq_str} > {fq_} 2>> {log}''', shell=True)

elif reads_format == "sra":
    r1 = output_fq_list[0]
    r2 = output_fq_list[1]
 
    if reads_num == 1:
        fq = input_fq_list[0]
        sra_file = os.path.basename(fq)
        sp.run(f'''rm -rf {output_dir}/{sra_file}* 2>> {log}''', shell=True)
        sp.run(f'''rm -rf {output_dir}.{sra_file}.temp 2>> {log}''', shell=True)

        sp.run(
            f''' 
            fasterq-dump \
            --threads {threads} \
            --split-3 \
            --temp {output_dir}.{sra_file}.temp \
            --outdir {output_dir} \
            {fq} 2>>{log}
            ''', shell=True)

        sp.run(f'''rm -rf {output_dir}.{sra_file}.temp 2>> {log}''', shell=True)
        sp.run(f'''pigz -f -p {threads} {output_dir}/{sra_file}_1.fastq 2>> {log}''', shell=True)
        sp.run(f'''pigz -f -p {threads} {output_dir}/{sra_file}_2.fastq 2>> {log}''', shell=True)
        sp.run(f'''rm -rf {output_dir}/{sra_file}._*.fastq 2>> {log}''', shell=True)

        sp.run(f'''mv {output_dir}/{sra_file}_1.fastq.gz {r1} 2>> {log}''', shell=True)
        sp.run(f'''mv {output_dir}/{sra_file}_2.fastq.gz {r2} 2>> {log}''', shell=True)

    else:
        r1_list = []
        r2_list = []
        for sra in input:
            sra_file = os.path.basename(sra)
            sra_1 = os.path.join(output_dir, sra_file + "_1.fastq.gz")
            sra_2 = os.path.join(output_dir, sra_file + "_2.fastq.gz")
            r1_list.append(sra_1)
            r2_list.append(sra_2)

            sp.run(f'''rm -rf {output_dir}/{sra_file}* 2>> {log}''', shell=True)
            sp.run(f'''rm -rf {output_dir}.{sra_file}.temp 2>> {log}''', shell=True)

            sp.run(
                f'''
                fasterq-dump \
                --threads {threads} \
                --split-3 \
                --temp {output_dir}.{sra_file}.temp \
                --outdir {output_dir} \
                {sra} 2>> {log}
                ''', shell=True)

            sp.run(f'''rm -rf {output_dir}.{sra_file}.temp 2>> {log}''', shell=True)
            sp.run(f'''pigz -f -p {threads} {output_dir}/{sra_file}_1.fastq 2>> {log}''', shell=True)
            sp.run(f'''pigz -f -p {threads} {output_dir}/{sra_file}_2.fastq 2>> {log}''', shell=True)
            sp.run(f'''rm -rf {output_dir}/{sra_file}._*.fastq 2>> {log}''', shell=True)

        r1_str = " ".join(r1_list)
        r2_str = " ".join(r2_list)
        sp.run(f'''cat {r1_str} > {r1} 2>> {log}''', shell=True)
        sp.run(f'''cat {r2_str} > {r2} 2>> {log}''', shell=True)
        sp.run(f'''rm -rf {r1_str} 2>> {log}''', shell=True)
        sp.run(f'''rm -rf {r2_str} 2>> {log}''', shell=True)