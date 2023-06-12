#!/usr/bin/env python

import os
import filecmp
import json
from executor import execute
from pprint import pprint


sample_id = str(snakemake.params.sample_id)
input_files = snakemake.params.input_files
headers = snakemake.params.headers
check_paired = snakemake.params.check_paired

threads = int(snakemake.threads)
log = str(snakemake.log)


pprint(input_files)

output = str(snakemake.output)
outdir = os.path.dirname(output)

outdir_pe = os.path.join(outdir, "pe")
outdir_se = os.path.join(outdir, "se")
outdir_long = os.path.join(outdir, "long")

outdir_pe_temp = os.path.join(outdir, "pe_temp")
outdir_se_temp = os.path.join(outdir, "se_temp")
outdir_long_temp = os.path.join(outdir, "long_temp")

is_fastq = True
input_tags = input_files.keys()
for k in input_tags:
    if k in headers["SRA"].values():
        is_fastq = False

# final file
r1 = os.path.join(outdir_pe, f"{sample_id}.raw.pe.1.fq.gz")
r2 = os.path.join(outdir_pe, f"{sample_id}.raw.pe.2.fq.gz")
rs = os.path.join(outdir_se, f"{sample_id}.raw.se.fq.gz")
rl = os.path.join(outdir_long, f"{sample_id}.raw.long.fq.gz")

# temp file
execute(f'''rm -rf {outdir}''')
execute(f'''mkdir -p {outdir}''')
execute(f'''rm -rf {log}''')

samples_dict ={}

if is_fastq:
    if headers["FQ"]["PE_FORWARD"] in input_tags:
        samples_dict["PE_FORWARD"] = r1
        samples_dict["PE_REVERSE"] = r2

        execute(f'''mkdir -p {outdir_pe}''')
        execute(f'''mkdir -p {outdir_pe_temp}''')

        r1_temp = os.path.join(outdir_pe_temp, f"{sample_id}.raw.temp.1.gz")
        r2_temp = os.path.join(outdir_pe_temp, f"{sample_id}.raw.temp.2.gz")
        r1_temp_p = os.path.join(outdir_pe_temp, f"{sample_id}.raw.temp.try_pe.1.gz")
        r2_temp_p = os.path.join(outdir_pe_temp, f"{sample_id}.raw.temp.try_pe.2.gz")

        forward_reads = input_files[headers["FQ"]["PE_FORWARD"]]
        reverse_reads = input_files[headers["FQ"]["PE_REVERSE"]]
        if len(forward_reads) == 1:
            fq1 = os.path.realpath(forward_reads[0])
            fq2 = os.path.realpath(reverse_reads[0])
            execute(f'''ln -s {fq1} {r1_temp} 2>> {log}''')
            execute(f'''ln -s {fq2} {r2_temp} 2>> {log}''')
        elif len(forward_reads) > 1:
            fq1 = " ".join(forward_reads)
            fq2 = " ".join(reverse_reads)
            execute(f'''cat {fq1} > {r1_temp} 2> {log}''')
            execute(f'''cat {fq2} > {r2_temp} 2>> {log}''')

        if check_paired:
            print("checking paired")
            id1 = f"{outdir_pe_temp}/id.list.1"
            id2 = f"{outdir_pe_temp}/id.list.2"
            execute(f'''seqkit seq -ni {r1_temp} | sed 's#/1$##g' > {id1} 2>> {log}''')
            execute(f'''seqkit seq -ni {r2_temp} | sed 's#/2$##g' > {id2} 2>> {log}''')

            if filecmp.cmp(id1, id2):
                execute(f'''mv {r1_temp} {r1} 2>> {log}''')
                execute(f'''mv {r2_temp} {r2} 2>> {log}''')
                execute(f'''rm -rf {id1}''')
                execute(f'''rm -rf {id2}''')
            else:
                idp = f"{outdir_pe_temp}/id.list.paired"
                execute(f'''mkdir -p {outdir_pe_temp}/temp0''')
                execute(f'''sort -T {outdir_pe_temp}/temp0 {id1} {id2} | uniq -c | awk '$1==2{{print $2}}' > {idp} 2>> {log}''')
                oneline = execute(f'''zcat {r1_temp} | head -1''', capture=True)
                if "/1" in oneline:
                    cmd1 = f'''awk '{{print $0 "/1"}}' {idp} | seqkit grep -f - {r1_temp} -o {r1_temp_p} 2>> {log}'''
                    cmd2 = f'''awk '{{print $0 "/2"}}' {idp} | seqkit grep -f - {r2_temp} -o {r2_temp_p} 2>> {log}'''
                    execute(cmd1)
                    execute(cmd2)
                else:
                    execute(f'''seqkit grep -f {idp} {r1_temp} -o {r1_temp_p} 2>> {log}''')
                    execute(f'''seqkit grep -f {idp} {r2_temp} -o {r2_temp_p} 2>> {log}''')

                execute(f'''rm -rf {r1_temp} 2>> {log}''')
                execute(f'''rm -rf {r2_temp} 2>> {log}''')
                execute(f'''rm -rf {id1} 2>> {log}''')
                execute(f'''rm -rf {id2} 2>> {log}''')
                execute(f'''rm -rf {idp} 2>> {log}''')

                ## more check
                execute(f'''mkdir -p {outdir_pe_temp}/temp1''')
                execute(f'''mkdir -p {outdir_pe_temp}/temp2''')
                execute(
                    f'''
                    pigz -b 102400 -dc {r1_temp_p} | \
                    paste - - - - | \
                    LC_ALL=C sort --parallel 4 -n -T {outdir_pe_temp}/temp1 -S 8G | \
                    tr '\t' '\n' | \
                    pigz -b 102400 > {r1} \
                    2>> {log}''')
                execute(
                    f'''
                    pigz -b 102400 -dc {r2_temp_p} | \
                    paste - - - - | \
                    LC_ALL=C sort --parallel 4 -n -T {outdir_pe_temp}/temp2 -S 8G | \
                    tr '\t' '\n' | \
                    pigz -b 102400 > {r2} \
                    2>> {log}''')

                execute(f'''rm -rf {r1_temp_p} 2>> {log}''')
                execute(f'''rm -rf {r2_temp_p} 2>> {log}''')
        else:
            execute(f'''mv {r1_temp} {r1} 2>>{log}''')
            execute(f'''mv {r2_temp} {r2} 2>>{log}''')

        execute(f'''rm -rf {outdir_pe_temp}''')

    elif headers["FQ"]["INTERLEAVED"] in input_tags:
        samples_dict["PE_FORWARD"] = r1
        samples_dict["PE_REVERSE"] = r2

        execute(f'''mkdir -p {outdir_pe}''')

        fq = " ".join(input_files[headers["FQ"]["INTERLEAVED"]])
        execute(
            f'''
            cat {fq} | \
            tee >(seqtk seq -1 - | pigz -cf -p {threads} > {r1}) | \
            seqtk seq -2 - | pigz -cf -p {threads} > {r2} 2>> {log}
            ''')

    if headers["FQ"]["SE"] in input_tags:
        samples_dict["SE"] = rs

        execute(f'''mkdir -p {outdir_se}''')

        single_reads = input_files[headers["FQ"]["SE"]]
        if len(single_reads) == 1:
            fq = os.path.realpath(single_reads[0])
            execute(f'''ln -s {fq} {rs} 2>> {log}''')
        elif len(single_reads) > 1:
            fq = " ".join(single_reads)
            execute(f'''cat {fq} > {rs} 2>> {log}''')

    if headers["FQ"]["LONG"] in input_tags:
        samples_dict["LONG"] = rl

        execute(f'''mkdir -p {outdir_long}''')

        long_reads = input_files[headers["FQ"]["LONG"]]
        if len(long_reads) == 1:
            fq = os.path.realpath(long_reads[0])
            execute(f'''ln -s {fq} {rl} 2>> {log}''')
        elif len(long_reads) > 1:
            fq = " ".join(long_reads)
            execute(f'''cat {fq} > {rl} 2>> {log}''')

else:
    if headers["SRA"]["PE"] in input_tags:
        samples_dict["PE_FORWARD"] = r1
        samples_dict["PE_REVERSE"] = r2

        execute(f'''mkdir -p {outdir_pe}''')
        execute(f'''mkdir -p {outdir_pe_temp}''')

        sra_pe = input_files[headers["SRA"]["PE"]]

        if len(sra_pe) == 1:
            sra = sra_pe[0]
            sra_name = os.path.basename(sra)
            execute(f'''rm -rf {outdir_pe_temp}/{sra_name}* 2>> {log}''')

            execute(
                f'''
                fasterq-dump \
                --threads {threads} \
                --split-3 \
                --temp {outdir_pe_temp}/{sra_name}.temp \
                --outdir {outdir_pe_temp} \
                {sra} 2>>{log}
                ''')

            execute(f'''rm -rf {outdir_pe_temp}/{sra_name}.temp* 2>> {log}''')
            execute(f'''pigz -f -p {threads} {outdir_pe_temp}/{sra_name}_1.fastq 2>> {log}''')
            execute(f'''pigz -f -p {threads} {outdir_pe_temp}/{sra_name}_2.fastq 2>> {log}''')
            execute(f'''rm -rf {outdir_pe_temp}/{sra_name}_*.fastq 2>> {log}''')
            execute(f'''mv {outdir_pe_temp}/{sra_name}_1.fastq.gz {r1} 2>> {log}''')
            execute(f'''mv {outdir_pe_temp}/{sra_name}_2.fastq.gz {r2} 2>> {log}''')

        else:
            r1_list = []
            r2_list = []
            for sra in sra_pe:
                sra_name = os.path.basename(sra)
                sra_1 = os.path.join(outdir_pe_temp, sra_name + "_1.fastq.gz")
                sra_2 = os.path.join(outdir_pe_temp, sra_name + "_2.fastq.gz")
                r1_list.append(sra_1)
                r2_list.append(sra_2)

                execute(f'''rm -rf {outdir_pe_temp}/{sra_name}* 2>> {log}''')

                execute(
                    f'''
                    fasterq-dump \
                    --threads {threads} \
                    --split-3 \
                    --temp {outdir_pe_temp}/{sra_name}.temp \
                    --outdir {outdir_pe_temp} \
                    {sra} 2>> {log}
                    ''')

                execute(f'''rm -rf {outdir_pe_temp}/{sra_name}.temp* 2>> {log}''')
                execute(f'''pigz -f -p {threads} {outdir_pe_temp}/{sra_name}_1.fastq 2>> {log}''')
                execute(f'''pigz -f -p {threads} {outdir_pe_temp}/{sra_name}_2.fastq 2>> {log}''')
                execute(f'''rm -rf {outdir_pe_temp}/{sra_name}._*.fastq 2>> {log}''')

            r1_str = " ".join(r1_list)
            r2_str = " ".join(r2_list)
            execute(f'''cat {r1_str} > {r1} 2>> {log}''')
            execute(f'''cat {r2_str} > {r2} 2>> {log}''')
            execute(f'''rm -rf {r1_str} 2>> {log}''')
            execute(f'''rm -rf {r2_str} 2>> {log}''')

        execute(f'''rm -rf {outdir_pe_temp}''')


    if headers["SRA"]["SE"] in input_tags:
        samples_dict["SE"] = rs

        execute(f'''mkdir -p {outdir_se}''')
        execute(f'''mkdir -p {outdir_se_temp}''')

        sra_se = input_files[headers["SRA"]["SE"]]

        if len(sra_se) == 1:
            sra = sra_se[0]
            sra_name = os.path.basename(sra)
            execute(f'''rm -rf {outdir_se_temp}/{sra_name}* 2>> {log}''')

            execute(
                f'''
                fasterq-dump \
                --threads {threads} \
                --temp {outdir_se_temp}/{sra_name}.temp \
                --outdir {outdir_se_temp} \
                {sra} 2>>{log}
                ''')

            execute(f'''rm -rf {outdir_se_temp}/{sra_name}.temp* 2>> {log}''')
            execute(f'''pigz -f -p {threads} {outdir_se_temp}/{sra_name}.fastq 2>> {log}''')
            execute(f'''rm -rf {outdir_se_temp}/{sra_name}*fastq 2>> {log}''')
            execute(f'''mv {outdir_se_temp}/{sra_name}.fastq.gz {rs} 2>> {log}''')

        else:
            rs_list = []
            for sra in sra_se:
                sra_name = os.path.basename(sra)
                r = os.path.join(outdir_se_temp, sra_name + ".fastq.gz")
                rs_list.append(r)

                execute(f'''rm -rf {outdir_se_temp}/{sra_name}* 2>> {log}''')

                execute(
                    f'''
                    fasterq-dump \
                    --threads {threads} \
                    --temp {outdir_se_temp}/{sra_name}.temp \
                    --outdir {outdir_se_temp} \
                    {sra} 2>> {log}
                    ''')

                execute(f'''rm -rf {outdir_se_temp}/{sra_name}.temp* 2>> {log}''')
                execute(f'''pigz -f -p {threads} {outdir_se_temp}/{sra_name}.fastq 2>> {log}''')
                execute(f'''rm -rf {outdir_se_temp}/{sra_name}*fastq 2>> {log}''')

            rs_str = " ".join(rs_list)
            execute(f'''cat {rs_str} > {rs} 2>> {log}''')
            execute(f'''rm -rf {rs_str} 2>> {log}''')

        execute(f'''rm -rf {outdir_se_temp}''')


    if headers["SRA"]["LONG"] in input_tags:
        samples_dict["LONG"] = rl

        execute(f'''mkdir -p {outdir_long}''')
        execute(f'''mkdir -p {outdir_long_temp}''')

        sra_l = input_files[headers["SRA"]["LONG"]]

        if len(sra_l) == 1:
            sra = sra_l[0]
            sra_name = os.path.basename(sra)
            execute(f'''rm -rf {outdir_long_temp}/{sra_name}* 2>> {log}''')

            execute(
                f'''
                fasterq-dump \
                --threads {threads} \
                --temp {outdir_long_temp}/{sra_name}.temp \
                --outdir {outdir_long_temp} \
                {sra} 2>>{log}
                ''')

            execute(f'''rm -rf {outdir_long_temp}/{sra_name}.temp* 2>> {log}''')
            execute(f'''pigz -f -p {threads} {outdir_long_temp}/{sra_name}.fastq 2>> {log}''')
            execute(f'''rm -rf {outdir_long_temp}/{sra_name}*fastq 2>> {log}''')
            execute(f'''mv {outdir_long_temp}/{sra_name}.fastq.gz {rl} 2>> {log}''')

        else:
            rl_list = []
            for sra in sra_l:
                sra_name = os.path.basename(sra)
                r = os.path.join(outdir_long_temp, sra_name + ".fastq.gz")
                rl_list.append(r)

                execute(f'''rm -rf {outdir_long_temp}/{sra_name}* 2>> {log}''')

                execute(
                    f'''
                    fasterq-dump \
                    --threads {threads} \
                    --temp {outdir_long_temp}/{sra_name}.temp \
                    --outdir {outdir_long_temp} \
                    {sra} 2>> {log}
                    ''')

                execute(f'''rm -rf {outdir_long_temp}/{sra_name}.temp* 2>> {log}''')
                execute(f'''pigz -f -p {threads} {outdir_long_temp}/{sra_name}.fastq 2>> {log}''')
                execute(f'''rm -rf {outdir_long_temp}/{sra_name}*fastq 2>> {log}''')

            rl_str = " ".join(rl_list)
            execute(f'''cat {rl_str} > {rl} 2>> {log}''')
            execute(f'''rm -rf {rl_str} 2>> {log}''')

        execute(f'''rm -rf {outdir_long_temp}''')


samples_obct = json.dumps(samples_dict, indent=2)
with open(output, 'wt') as oh:
    oh.write(samples_obct)
