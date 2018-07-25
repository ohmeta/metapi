rule coverage_metabat2:
    input:
        bam = os.path.join(config["results"]["alignment"], "{sample}.sorted.bam")
    output:
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.metabat2.depth.txt")
    params:
        depth_dir = config["results"]["binning"]["depth"]
    shell:
        '''
        mkdir -p {params.depth_dir}
        docker run metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        '''

rule coverage_maxbin2:
    input:
        bam = os.path.join(config["results"]["alignment"], "{sample}.sorted.bam")
    output:
        depth_bb = os.path.join(config["results"]["binning"]["depth"], "{sample}.bbmap.depth.txt"),
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.maxbin2.depth.txt")
    shell:
        '''
        pileup.sh in={input.bam} out={output.depth_bb}
        awk '{print $1 "\t" $5}' {output.depth_bb} | grep -v '^#' > {output.depth}
        '''

rule binning_metabat2:
    input:
        asmfa = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa"),
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.metabat2.depth.txt")
    output:
        default = os.path.join(config["logs"]["binning"]["metabat2"], "{sample}.metabat2.done")
    log:
        os.path.join(config["logs"]["binning"]["metabat2"], "{sample}.metabat2.log")
    params:
        bins_dir = os.path.join(config["results"]["binning"]["bins"], "{sample}.metabat2_out" ),
        bin_prefix = os.path.join(config["results"]["binning"]["bins"], "{sample}.metabat2_out/{sample}.bin"),
        min_contig = config["params"]["binning"]["metabat2"]["min_contig"],
        seed = config["params"]["binning"]["metabat2"]["seed"]
    shell:
        '''
        mkdir -p {params.bins_dir}
        docker run metabat/metabat:latest metabat2 -i {input.asmfa} -a {input.depth} -o {params.bin_prefix} -m {params.min_contig} --seed {params.seed} -v > {log}
        touch {output.default}
        '''

rule binning_maxbin2:
    input:
        asmfa = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa"),
        depth = os.path.join(config["results"]["binning"]["depth"], "{sample}.maxbin2.depth.txt")
    output:
        os.path.join(config["results"]["binning"]["bins"], "{sample}.maxbin2_out/{sample}.bin.summary")
    params:
        bins_dir = os.path.join(config["results"]["binning"]["bins"], "{sample}.maxbin2_out"),
        bin_prefix = os.path.join(config["results"]["binning"]["bins"], "{sample}.maxbin2_out/{sample}.bin")
    threads:
        config['params']['binning']['maxbin2']['threads']
    shell:
        '''
        mkdir -p {params.bins_dir}
        run_MaxBin.pl -thread {threads} -contig {input.asmfa} -out {params.bin_prefix} -abund {intput.depth}
        '''
