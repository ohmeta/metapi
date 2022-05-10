def rmhost_input(wildcards, have_single=False):
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming", have_single)
    else:
        return get_reads(wildcards, "raw", have_single)


def trimming_stats_input(wildcards, have_single=False):
    if TRIMMING_DO and config["params"]["qcreport"]["do"]:
        return expand(os.path.join(config["output"]["trimming"],
                                   "report/stats/{sample}_trimming_stats.tsv"),
                      sample=wildcards.sample)
    else:
        return []
 

if config["params"]["rmhost"]["soap"]["do"]:
    rule rmhost_soap_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["soap"]["index_prefix"],
                   suffix=["amb", "ann", "bwt", "fmv", "hot", "lkt", "pac",
                           "rev.bwt", "rev.fmv", "rev.lkt", "rev.pac", "sa", ".sai"])
        log:
            os.path.join(config["output"]["rmhost"],
                         "logs/build_host_index_for_soap.log")
        params:
            host_fasta = os.path.basename(config["params"]["rmhost"]["host_fasta"]),
            index_dir = os.path.dirname(config["params"]["rmhost"]["soap"]["index_prefix"])
        shell:
            '''
            mkdir -p {params.index_dir}
            pushd {params.index_dir} && \
            ln -s {input} && \
            2bwt-builder {params.host_fasta} 2> {log} && popd
            '''


    rule rmhost_soap:
        input:
            lambda wildcards: trimming_stats_input(wildcards),
            reads = lambda wildcards: rmhost_input(wildcards),
            index = expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["soap"]["index_prefix"],
                   suffix=["amb", "ann", "bwt", "fmv", "hot", "lkt", "pac",
                           "rev.bwt", "rev.fmv", "rev.lkt", "rev.pac", "sa", "sai"])
        output:
            reads = expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "") \
                           if config["params"]["rmhost"]["save_reads"] else \
                              temp(expand(os.path.join(
                                  config["output"]["rmhost"],
                                  "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                          read=[".1", ".2"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.soap.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/soap/{sample}.soap.txt")
        priority:
            10
        params:
            index_prefix = config["params"]["rmhost"]["soap"]["index_prefix"],
            match_model = config["params"]["rmhost"]["soap"]["match_model"],
            align_seed = config["params"]["rmhost"]["soap"]["align_seed"],
            report_repeat_hits = config["params"]["rmhost"]["soap"]["report_repeat_hits"],
            max_mismatch_num = config["params"]["rmhost"]["soap"]["max_mismatch_num"],
            identity = config["params"]["rmhost"]["soap"]["identity"]
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            outdir = os.path.dirname(output.reads[0])
            shell(f'''mkdir -p {outdir}''')

            if IS_PE:
                shell(
                    f'''
                    soap2.22 \
                    -a {input.reads[0]} -b {input.reads[1]} \
                    -D {params.index_prefix} \
                    -M {params.match_model} \
                    -l {params.align_seed} \
                    -r {params.report_repeat_hits} \
                    -v {params.max_mismatch_num} \
                    -c {params.identity} \
                    -p {threads} \
                    -S \
                    -o {outdir}/soap.pe \
                    -2 {outdir}/soap.se \
                    -u {outdir}/unmapped.fq 2> {log}

                    seqtk dropse {outdir}/unmapped.fq |
                    tee >(seqtk seq -1 - | pigz > {output.reads[0]}) |
                    seqtk seq -2 - | pigz > {output.reads[1]}

                    rm -rf {outdir}/soap.pe
                    rm -rf {outdir}/soap.se
                    rm -rf {outdir}/unmapped.fq
                    ''')
            else:
                shell(
                    f'''
                    soap2.22 \
                    -a {input.reads[0]} \
                    -D {params.index_prefix} \
                    -M {params.match_model} \
                    -l {params.align_seed} \
                    -r {params.report_repeat_hits} \
                    -v {params.max_mismatch_num} \
                    -c {params.identity} \
                    -p {threads} \
                    -S \
                    -o {outdir}/soap.se \
                    -u {outdir}/unmapped.fq 2> {log}

                    rm -rf {outdir}/soap.se
                    pigz {outdir}/unmapped.fq
                    mv {outdir}/unmapped.fq.gz {output.reads[0]}
                    ''')


    rule rmhost_soap_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_soap_all:
        input:


BWA_INDEX_SUFFIX = ["0123", "amb", "ann", "bwt.2bit.64", "pac"] if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" \
else ["amb", "ann", "bwt", "pac", "sa"]

if config["params"]["rmhost"]["bwa"]["do"]:
    rule rmhost_bwa_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                   suffix=BWA_INDEX_SUFFIX)
        benchmark:
            os.path.join(config["output"]["rmhost"],
                         "benchmark/bwa.index.benchmark.txt")
        log:
            os.path.join(config["output"]["rmhost"],
                         "logs/build_host_index_for_bwa.log")
        params:
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bwa = "bwa-mem2" if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" else "bwa"
        shell:
            '''
            {params.bwa} index -p {params.index_prefix} {input} 2> {log}
            '''


    rule rmhost_bwa:
        input:
            lambda wildcards: trimming_stats_input(wildcards),
            reads = lambda wildcards: rmhost_input(wildcards),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bwa"]["index_prefix"],
                           suffix=BWA_INDEX_SUFFIX)
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/flagstat/{sample}.align2host.flagstat"),
            reads = expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "") \
                           if config["params"]["rmhost"]["save_reads"] else \
                              temp(expand(os.path.join(
                                  config["output"]["rmhost"],
                                  "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                          read=[".1", ".2"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.bwa.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/bwa/{sample}.bwa.txt")
        priority:
            10
        params:
            compression = config["params"]["rmhost"]["compression"],
            bwa = "bwa-mem2" if config["params"]["rmhost"]["bwa"]["algorithms"] == "mem2" else "bwa",
            minimum_seed_length = config["params"]["rmhost"]["bwa"]["minimum_seed_length"],
            index_prefix = config["params"]["rmhost"]["bwa"]["index_prefix"],
            bam = os.path.join(config["output"]["rmhost"],
                               "bam/{sample}/{sample}.align2host.sorted.bam")
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["save_bam"]:
                    shell("mkdir -p %s" % os.path.dirname(params.bam))
                    shell(
                        '''
                        {params.bwa} mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} {input.reads[1]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 12 -F 256 \
                              -1 {output.reads[0]} \
                              -2 {output.reads[1]} -) | \
                        samtools sort \
                        -m 3G \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} - \
                        2> {log}
                        ''')
                else:
                    shell(
                        '''
                        {params.bwa} mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} {input.reads[1]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 {output.reads[0]} \
                        -2 {output.reads[1]} - \
                        2> {log}
                        ''')
            else:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        {params.bwa} mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 4 -F 256 - | \
                              pigz -c -p {threads} \
                              > {output.reads[0]}) | \
                        samtools sort \
                        -m 3G \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} - \
                        2> {log}
                        ''')
                else:
                    shell(
                        '''
                        {params.bwa} mem \
                        -k {params.minimum_seed_length} \
                        -t {threads} \
                        {params.index_prefix} \
                        {input.reads[0]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -c -p {threads} \
                        > {output.reads[0]} \
                        2> {log}
                        ''')


    rule rmhost_bwa_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_bwa_all:
        input:


if config["params"]["rmhost"]["bowtie2"]["do"]:
    rule rmhost_bowtie2_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            expand("{prefix}.{suffix}",
                   prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                   suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        log:
            os.path.join(config["output"]["rmhost"], "logs/build_host_index_for_bowtie2.log")
        params:
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"]
        shell:
            '''
            bowtie2-build {input} {params.index_prefix} 2> {log}
            '''


    rule rmhost_bowtie2:
        input:
            lambda wildcards: trimming_stats_input(wildcards),
            reads = lambda wildcards: rmhost_input(wildcards),
            index = expand("{prefix}.{suffix}",
                           prefix=config["params"]["rmhost"]["bowtie2"]["index_prefix"],
                           suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/flagstat/{sample}.align2host.flagstat"),
            reads = expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "") \
                           if config["params"]["rmhost"]["save_reads"] else \
                              temp(expand(os.path.join(
                                  config["output"]["rmhost"],
                                  "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                          read=[".1", ".2"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.bowtie2.log")
        benchmark:
            os.path.join(config["output"]["rmhost"], "benchmark/bowtie2/{sample}.bowtie2.txt")
        priority:
            10
        params:
            presets = config["params"]["rmhost"]["bowtie2"]["presets"],
            compression = config["params"]["rmhost"]["compression"],
            index_prefix = config["params"]["rmhost"]["bowtie2"]["index_prefix"],
            bam = os.path.join(config["output"]["rmhost"],
                               "bam/{sample}/{sample}.align2host.sorted.bam")
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        bowtie2 \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -1 {input.reads[0]} \
                        -2 {input.reads[1]} \
                        {params.presets} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 12 -F 256 \
                              -1 {output.reads[0]} \
                              -2 {output.reads[1]} -) | \
                        samtools sort \
                        -m 3G \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} -
                        ''')
                else:
                    shell(
                        '''
                        bowtie2 \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -1 {input.reads[0]} \
                        -2 {input.reads[1]} \
                        {params.presets} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 {output.reads[0]} \
                        -2 {output.reads[1]} -
                        ''')
            else:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        bowtie2 \
                        {params.presets} \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -U {input.reads[0]} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 4 -F 256 - | \
                              pigz -c -p {threads} \
                              > {output.reads[0]}) | \
                        samtools sort \
                        -m 3G \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} -
                        ''')
                else:
                    shell(
                        '''
                        bowtie2 \
                        {params.presets} \
                        --threads {threads} \
                        -x {params.index_prefix} \
                        -U {input.reads[0]} \
                        2> {log} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -c -p {threads} \
                        > {output.reads[0]}
                        ''')


    rule rmhost_bowtie2_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_bowtie2_all:
        input:


if config["params"]["rmhost"]["minimap2"]["do"]:
    rule rmhost_minimap2_index:
        input:
            config["params"]["rmhost"]["host_fasta"]
        output:
            config["params"]["rmhost"]["minimap2"]["index"]
        params:
            split_size = config["params"]["rmhost"]["minimap2"]["split_size"]
        shell:
            '''
            minimap2 -I {params.split_size} -d {output} {intput}
            ''' 


    rule rmhost_minimap2:
        input:
            lambda wildcards: trimming_stats_input(wildcards),
            reads = lambda wildcards: rmhost_input(wildcards),
            index = config["params"]["rmhost"]["minimap2"]["index"] 
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/flagstat/{sample}.align2host.flagstat"),
            reads = expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "") \
                           if config["params"]["rmhost"]["save_reads"] else \
                              temp(expand(os.path.join(
                                  config["output"]["rmhost"],
                                  "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                          read=[".1", ".2"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.minimap2.log")
        benchmark:
            os.path.join(config["output"]["rmhost"],
                         "benchmark/minimap2/{sample}.minimap2.txt")
        priority:
            10
        params:
            preset = config["params"]["rmhost"]["minimap2"]["preset"],
            compression = config["params"]["rmhost"]["compression"],
            bam = os.path.join(config["output"]["rmhost"],
                               "bam/{sample}/{sample}.align2host.sorted.bam")
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            if IS_PE:
                if config["params"]["rmhost"]["save_bam"]:
                    shell("mkdir -p %s" % os.path.dirname(params.bam))
                    shell(
                        '''
                        minimap2 \
                        -t {threads} \
                        -ax {params.preset} \
                        {input.index} \
                        {input.reads[0]} {input.reads[1]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 12 -F 256 \
                              -1 {output.reads[0]} \
                              -2 {output.reads[1]} -) | \
                        samtools sort \
                        -m 3G \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} - \
                        2> {log}
                        ''')
                else:
                    shell(
                        '''
                        minimap2 \
                        -t {threads} \
                        -ax {params.preset} \
                        {input.index} \
                        {input.reads[0]} {input.reads[1]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 12 -F 256 \
                        -1 {output.reads[0]} \
                        -2 {output.reads[1]} - \
                        2> {log}
                        ''')
            else:
                if config["params"]["rmhost"]["save_bam"]:
                    shell('''mkdir -p %s''' % os.path.dirname(params.bam))
                    shell(
                        '''
                        minimap2 \
                        -t {threads} \
                        {input.index} \
                        {input.reads[0]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        tee >(samtools fastq \
                              -@{threads} \
                              -c {params.compression} \
                              -N -f 4 -F 256 - | \
                              pigz -c -p {threads} \
                              > {output.reads[0]}) | \
                        samtools sort \
                        -m 3G \
                        -@{threads} \
                        -T {params.bam} \
                        -O BAM -o {params.bam} - \
                        2> {log}
                        ''')
                else:
                    shell(
                        '''
                        minimap2 \
                        -t {threads} \
                        {input.index} \
                        {input.reads[0]} | \
                        tee >(samtools flagstat \
                              -@{threads} - \
                              > {output.flagstat}) | \
                        samtools fastq \
                        -@{threads} \
                        -c {params.compression} \
                        -N -f 4 -F 256 - | \
                        pigz -c -p {threads} \
                        > {output.reads[0]} \
                        2> {log}
                        ''')


    rule rmhost_minimap2_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_minimap2_all:
        input:


if config["params"]["rmhost"]["kraken2"]["do"]:
    rule rmhost_kraken2:
        input:
            lambda wildcards: trimming_stats_input(wildcards),
            reads = lambda wildcards: rmhost_input(wildcards),
            database = config["params"]["rmhost"]["kraken2"]["database"] 
        output:
            table = temp(os.path.join(config["output"]["rmhost"],
                                      "short_reads/{sample}/{sample}.kraken2.table")),
            report = os.path.join(config["output"]["rmhost"],
                                  "short_reads/{sample}/{sample}.kraken2.report.gz"),
            reads = expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "") \
                           if config["params"]["rmhost"]["save_reads"] else \
                              temp(expand(os.path.join(
                                  config["output"]["rmhost"],
                                  "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                          read=[".1", ".2"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.kraken2.log")
        benchmark:
            os.path.join(config["output"]["rmhost"],
                         "benchmark/kraken2/{sample}.kraken2.txt")
        params:
            confidence = config["params"]["rmhost"]["kraken2"]["confidence"],
            min_base_quality = config["params"]["rmhost"]["kraken2"]["min_base_quality"],
            min_hit_groups = config["params"]["rmhost"]["kraken2"]["min_hit_groups"],
            host_taxid = config["params"]["rmhost"]["kraken2"]["host_taxid"]
        priority:
            10
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            import os
            report = os.path.splitext(output.report)[0]

            if IS_PE:
                shell(
                    '''
                    kraken2 \
                    --threads {threads} \
                    --db {input.database} \
                    --use-names \
                    --confidence {params.confidence} \
                    --minimum-base-quality {params.min_base_quality} \
                    --minimum-hit-groups {params.min_hit_groups} \
                    --output {output.table} \
                    --report {report} \
                    --gzip-compressed \
                    --paired \
                    {input.reads} \
                    >{log} 2>&1

                    pigz -p {threads} {report}

                    extract_kraken2_reads.py \
                    -k {output.table} \
                    --taxid {params.host_taxid} \
                    --noappend \
                    --exclude \
                    --fastq-output \
                    --gzip-output \
                    -s {input.reads[0]} \
                    -s2 {input.reads[1]} \
                    -o {output.reads[0]} \
                    -o2 {output.reads[1]} \
                    >>{log} 2>&1
                    ''')
            else:
                shell(
                    '''
                    kraken2 \
                    --threads {threads} \
                    --db {input.database} \
                    --use-names \
                    --confidence {params.confidence} \
                    --minimum-base-quality {params.min_base_quality} \
                    --minimum-hit-groups {params.min_hit_groups} \
                    --output {output.table} \
                    --report {report} \
                    --gzip-compressed \
                    {input.reads[0]} \
                    >{log} 2>&1

                    pigz -p {threads} {report}

                    extract_kraken2_reads.py \
                    -k {output.table} \
                    --taxid {params.host_taxid} \
                    --noappend \
                    --exclude \
                    --fastq-output \
                    --gzip-output \
                    -s {input.reads[0]} \
                    -o {output.reads[0]} \
                    >>{log} 2>&1
                    ''')


    rule rmhost_kraken2_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_kraken2_all:
        input:


if config["params"]["rmhost"]["kneaddata"]["do"]:
    rule rmhost_kneaddata:
        input:
            lambda wildcards: trimming_stats_input(wildcards),
            reads = lambda wildcards: rmhost_input(wildcards)
        output:
            reads = expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                           read=[".1", ".2"] if IS_PE else "") \
                           if config["params"]["rmhost"]["save_reads"] else \
                              temp(expand(os.path.join(
                                  config["output"]["rmhost"],
                                  "short_reads/{{sample}}/{{sample}}.rmhost{read}.fq.gz"),
                                          read=[".1", ".2"] if IS_PE else ""))
        log:
            os.path.join(config["output"]["rmhost"], "logs/{sample}.kneaddata.log")
        benchmark:
            os.path.join(config["output"]["rmhost"],
                         "benchmark/kneaddata/{sample}.kneaddata.txt")
        params:
            trf_options = "--run-trf" if config["params"]["rmhost"]["kneaddata"]["do_trf"] else "--bypass-trf",
            do_trimmomatic = config["params"]["rmhost"]["kneaddata"]["do_trimmomatic"],
            trimmomatic_options = config["params"]["rmhost"]["kneaddata"]["trimmomatic_options"],
            sequencer_source = config["params"]["rmhost"]["kneaddata"]["sequencer_source"],
            do_bowtie2 = config["params"]["rmhost"]["kneaddata"]["do_bowtie2"],
            bowtie2_options = config["params"]["rmhost"]["kneaddata"]["bowtie2_options"],
            decontaminate_pairs = config["params"]["rmhost"]["kneaddata"]["decontaminate_pairs"],
            bowtie2_database = config["params"]["rmhost"]["kneaddata"]["bowtie2_database"],
            do_bmtagger = config["params"]["rmhost"]["kneaddata"]["do_bmtagger"],
            output_dir = os.path.join(config["output"]["rmhost"], "short_reads/{sample}"),
            output_prefix = "{sample}.rmhost"
        priority:
            10
        threads:
            config["params"]["rmhost"]["threads"]
        run:
            import shutil
            import os
            trimmomatic_dir = os.path.dirname(os.path.realpath(shutil.which("trimmomatic")))

            shell('''rm -rf {params.output_dir}''')

            if params.do_bowtie2:
                if params.do_trimmomatic:
                    shell(
                        '''
                        kneaddata %s \
                        {params.trf_options} \
                        --output {params.output_dir} \
                        --output-prefix {params.output_prefix} \
                        --reference-db {params.bowtie2_database} \
                        --trimmomatic {trimmomatic_dir} \
                        --trimmomatic-options '{params.trimmomatic_options}' \
                        --sequencer-source {params.sequencer_source} \
                        --bowtie2-options '{params.bowtie2_options} ' \
                        --decontaminate-pairs {params.decontaminate_pairs} \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                        ''' % "-i {input.reads[0]} -i {input.reads[1]}" if IS_PE else "-i {input.reads}")

                else:
                    shell(
                        '''
                        kneaddata %s \
                        {params.trf_options} \
                        --bypass-trim \
                        --output {params.output_dir} \
                        --output-prefix {params.output_prefix} \
                        --reference-db {params.bowtie2_database} \
                        --bowtie2-options '{params.bowtie2_options} ' \
                        --decontaminate-pairs {params.decontaminate_pairs} \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                        ''' % "-i {input.reads[0]} -i {input.reads[1]}" if IS_PE else "-i {input.reads}")

            elif params.do_bmtagger:
                if params.do_trimmomatic:
                    shell(
                        '''
                        kneaddata %s \
                        {params.trf_options} \
                        --output {params.output_dir} \
                        --output-prefix {params.output_prefix} \
                        --trimmomatic {trimmomatic} \
                        --trimmomatic-options '{params.trimmomatic_options}' \
                        --sequencer-source {params.sequencer_source} \
                        --run-bmtagger \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                        ''' % "-i {input.reads[0]} -i {input.reads[1]}" if IS_PE else "-i {input.reads}")

                else:
                    shell(
                        '''
                        kneaddata %s \
                        {params.trf_options} \
                        --bypass-trim \
                        --output {params.output_dir} \
                        --output-prefix {params.output_prefix} \
                        --run-bmtagger \
                        --remove-intermediate-output \
                        --threads {threads} \
                        --reorder \
                        --log {log}
                        ''' % "-i {input.reads[0]} -i {input.reads[1]}" if IS_PE else "-i {input.reads}")
                 
            shell(
                '''
                pigz -p {threads} {params.output_dir}/* 
                ''')
            
            if IS_PE:
                shell(
                    '''
                    mv {params.output_dir}/{params.output_prefix}_paired_1.fastq.gz \
                    {output.reads[0]}
                    mv {params.output_dir}/{params.output_prefix}_paired_2.fastq.gz \
                    {output.reads[1]}
                    ''')
            else:
                shell(
                    '''
                    mv {params.output_dir}/{params.output_prefix}.fastq.gz \
                    {output.reads}
                    ''')


    rule rmhost_kneaddata_all:
        input:
            expand(os.path.join(
                config["output"]["rmhost"],
                "short_reads/{sample}/{sample}.rmhost{read}.fq.gz"),
                   sample=SAMPLES_ID_LIST,
                   read=[".1", ".2"] if IS_PE else "")

else:
    rule rmhost_kneaddata_all:
        input:


if RMHOST_DO \
and (not config["params"]["rmhost"]["soap"]["do"]) \
and (not config["params"]["rmhost"]["kraken2"]["do"]) \
and (not config["params"]["rmhost"]["kneaddata"]["do"]):
    rule rmhost_alignment_report:
        input:
            expand(
                os.path.join(config["output"]["rmhost"],
                             "report/flagstat/{sample}.align2host.flagstat"),
                sample=SAMPLES_ID_LIST)
        output:
            flagstat = os.path.join(config["output"]["rmhost"],
                                    "report/rmhost_align2host_stats.tsv")
        run:
            input_list = [str(i) for i in input]
            metapi.flagstats_summary(input_list, 2, output=output.flagstat)

else:
    rule rmhost_alignment_report:
        input:


if RMHOST_DO and config["params"]["qcreport"]["do"]:
    rule rmhost_report:
        input:
            lambda wildcards: get_reads(wildcards, "rmhost")
        output:
            os.path.join(config["output"]["rmhost"],
                         "report/stats/{sample}_rmhost_stats.tsv")
        log:
            os.path.join(config["output"]["rmhost"],
                         "logs/{sample}.seqkit.log")
        priority:
            25
        params:
            sample_id = "{sample}",
            fq_encoding = config["params"]["fq_encoding"]
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            shell(
                '''
                seqkit stats \
                --all \
                --basename \
                --tabular \
                --fq-encoding {params.fq_encoding} \
                --out-file {output} \
                --threads {threads} \
                {input} 2> {log}
                ''')

            if IS_PE:
                metapi.change(output[0], params.sample_id, "rmhost",
                              "pe", ["fq1", "fq2"])
            else:
                metapi.change(output[0], params.sample_id, "rmhost",
                              "se", ["fq1"])


    rule rmhost_report_merge:
        input:
            expand(
                os.path.join(config["output"]["rmhost"],
                             "report/stats/{sample}_rmhost_stats.tsv"),
                sample=SAMPLES_ID_LIST)
        output:
            os.path.join(config["output"]["qcreport"], "rmhost_stats.tsv")
        threads:
            config["params"]["qcreport"]["seqkit"]["threads"]
        run:
            metapi.merge(input, metapi.parse, threads, output=output[0])


    rule rmhost_report_all:
        input:
            os.path.join(config["output"]["qcreport"], "rmhost_stats.tsv")

else:
    rule rmhost_report_all:
        input:


rule rmhost_all:
    input:
        rules.rmhost_bwa_all.input,
        rules.rmhost_bowtie2_all.input,
        rules.rmhost_soap_all.input,
        rules.rmhost_minimap2_all.input,
        rules.rmhost_kraken2_all.input,
        rules.rmhost_kneaddata_all.input,
        rules.rmhost_alignment_report.input,
        rules.rmhost_report_all.input


localrules:
    rmhost_bwa_all,
    rmhost_soap_all,
    rmhost_bowtie2_all,
    rmhost_kraken2_all,
    rmhost_minimap2_all,
    rmhost_kneaddata_all,
    rmhost_report_all,
    rmhost_all