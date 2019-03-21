rule coverage_metabat2:
    input:
        # scaftigs_flagstat_summary = os.path.join(config["results"]["alignment"], "scaftigs_flagstat_summary.tsv"),
        bam = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam")
    output:
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.{assembler}.metabat2.depth.txt")
    params:
        depth_dir = directory(config["results"]["binning"]["depth"])
    shell:
        '''
        mkdir -p {params.depth_dir}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        '''

rule coverage_maxbin2:
    input:
        bam = os.path.join(config["results"]["alignment"], "{sample}.bwa_out/{sample}.{assembler}.sorted.bam")
    output:
        depth_bb = os.path.join(config["results"]["binning"]["depth"], "{sample}.{assembler}.bbmap.depth.txt"),
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.{assembler}.maxbin2.depth.txt")
    shell:
        '''
        pileup.sh in={input.bam} out={output.depth_bb}
        awk '{print $1 "\t" $5}' {output.depth_bb} | grep -v '^#' > {output.depth}
        '''

rule binning_metabat2:
    input:
        scaftigs = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz"),
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.{assembler}.metabat2.depth.txt")
    output:
        bins_dir = directory(os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.metabat2_out"))
    log:
        os.path.join(config["logs"]["binning"]["metabat2"], "{sample}.{assembler}.metabat2.log")
    params:
        bin_prefix = os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.metabat2_out/{sample}.{assembler}.bin"),
        min_contig = config["params"]["binning"]["metabat2"]["min_contig"],
        seed = config["params"]["binning"]["metabat2"]["seed"]
    shell:
        '''
        mkdir -p {output.bins_dir}
        metabat2 -i {input.scaftigs} -a {input.depth} -o {params.bin_prefix} -m {params.min_contig} --seed {params.seed} -v 2> {log}
        '''

rule binning_maxbin2:
    input:
        scaftigs = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz"),
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.{assembler}.maxbin2.depth.txt")
    output:
        os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.maxbin2_out/{sample}.{assembler}.bin.summary")
    params:
        bins_dir = directory(os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.maxbin2_out")),
        bin_prefix = os.path.join(config["results"]["binning"]["bins"], "{sample}.{assembler}.maxbin2_out/{sample}.{assembler}.bin")
    threads:
        config["params"]["binning"]["maxbin2"]["threads"]
    log:
        os.path.join(config["logs"]["binning"]["maxbin2"], "{sample}.{assembler}.maxbin2.log")
    shell:
        '''
        mkdir -p {params.bins_dir}
        run_MaxBin.pl -thread {threads} -contig {input.scaftigs} -out {params.bin_prefix} -abund {intput.depth} 2> {log}
        '''
'''
rule dastools:
    input:
    output:
    params:
    threads:
    log:
    shell:
'''

# ruleorder: summary_scaftigs_flagstat > coverage_metabat2