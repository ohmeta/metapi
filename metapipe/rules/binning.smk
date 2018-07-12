rule coverage:
    input:
        bam = os.path.join(config["results"]["alignment"], "{sample}.sorted.bam")
    output:
        depth = os.path.join(config["results"]["binning"], "coverage/{sample}.metabat2.depth.txt")
    params:
        coverage_dir = os.path.join(config["results"]["binning"], "coverage")
    shell:
        '''
        mkdir -p {params.coverage_dir}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        '''

rule metabat2:
    input:
        asmfa = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa"),
        depth = os.path.join(config["results"]["binning"], "coverage/{sample}.metabat2.depth.txt")
    output:
        default = os.path.join(config["logs"]["binning"], "{sample}.done")
    log:
        os.path.join(config["logs"]["binning"], "{sample}.metabat2.log")
    params:
        bins_dir = os.path.join(config["results"]["binning"], "bins"),
        bin_prefix = os.path.join(config["results"]["binning"], "bins/{sample}.metabat2_out/{sample}.bin"),
        min_contig = config["params"]["binning"]["metabat2"]["min_contig"],
        seed = config["params"]["binning"]["metabat2"]["seed"]
    shell:
        '''
        mkdir -p {params.bins_dir}
        metabat2 -i {input.asmfa} -a {input.depth} -o {params.bin_prefix} -m {params.min_contig} --seed {params.seed} -v > {log}
        touch {output.default}
        '''
