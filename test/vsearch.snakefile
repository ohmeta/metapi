'''vsearch v2.6.0_linux_x86_64, 15.6GB RAM, 8 cores
https://github.com/torognes/vsearch

For help, please enter: vsearch --help

For further details, please see the manual.

Example commands:

vsearch --allpairs_global FILENAME --id 0.5 --alnout FILENAME
vsearch --cluster_fast FILENAME --id 0.97 --centroids FILENAME
vsearch --cluster_size FILENAME --id 0.97 --centroids FILENAME
vsearch --cluster_smallmem FILENAME --usersort --id 0.97 --centroids FILENAME
vsearch --derep_fulllength FILENAME --output FILENAME
vsearch --derep_prefix FILENAME --output FILENAME
vsearch --fastq_chars FILENAME
vsearch --fastq_convert FILENAME --fastqout FILENAME --fastq_ascii 64
vsearch --fastq_eestats FILENAME --output FILENAME
vsearch --fastq_eestats2 FILENAME --output FILENAME
vsearch --fastq_mergepairs FILENAME --reverse FILENAME --fastqout FILENAME
vsearch --fastq_stats FILENAME --log FILENAME
vsearch --fastx_filter FILENAME --fastaout FILENAME --fastq_trunclen 100
vsearch --fastx_mask FILENAME --fastaout FILENAME
vsearch --fastx_revcomp FILENAME --fastqout FILENAME
vsearch --fastx_subsample FILENAME --fastaout FILENAME --sample_pct 1
vsearch --rereplicate FILENAME --output FILENAME
vsearch --search_exact FILENAME --db FILENAME --alnout FILENAME
vsearch --shuffle FILENAME --output FILENAME
vsearch --sortbylength FILENAME --output FILENAME
vsearch --sortbysize FILENAME --output FILENAME
vsearch --uchime_denovo FILENAME --nonchimeras FILENAME
vsearch --uchime_ref FILENAME --db FILENAME --nonchimeras FILENAME
vsearch --usearch_global FILENAME --db FILENAME --id 0.97 --alnout FILENAME'''

import os
r1 = {}
r2 = {}
for i in os.listdir("../data/00.raw/"):
    if i.endswith("1.fq.gz"):
        sn = i.rstrip(".1.fq.gz")
        r1[sn] = os.path.join("../data/00.raw", i)
    elif i.endswith("2.fq.gz"):
        sn = i.rstrip(".2.fq.gz")
        r2[sn] = os.path.join("../data/00.raw", i)

rule vsearch_test:
    input:
        reads_1 = lambda wildcards: r1[wildcards.sample]
        reads_2 = lambda wildcards: r2[wildcards.sample]
    output:

    shell:
        """
        vsearch --
        """
