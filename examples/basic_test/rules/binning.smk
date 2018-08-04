rule coverage:
    input:
        bam = os.path.join(config["results"]["alignment"], "{sample}_{unit}.sorted.bam")
    output:
        depth = os.path.join(config["results"]["binning"], "coverage/{sample}_{unit}.metabat2.depth.txt")
    params:
        coverage_dir = os.path.join(config["results"]["binning"], "coverage")
    shell:
        '''
        mkdir -p {params.coverage_dir}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        '''

rule metabat2:
    input:
        asmfa = os.path.join(config["results"]["assembly"], "{sample}_{unit}.megahit_out/{sample}_{unit}.contigs.fa"),
        depth = os.path.join(config["results"]["binning"], "coverage/{sample}_{unit}.metabat2.depth.txt")
    output:
        fa_dir = os.path.join(config["results"]["binning"], "bins/{sample}_{unit}.metabat2_out/"),
        default = os.path.join(config["logs"]["binning"], "{sample}_{unit}.done")
    log:
        os.path.join(config["logs"]["binning"], "{sample}_{unit}.metabat2.log")
    params:
        bins_dir = os.path.join(config["results"]["binning"], "bins"),
        bin_prefix = os.path.join(config["results"]["binning"], "bins/{sample}_{unit}.metabat2_out/{sample}_{unit}.bin"),
        min_contig = config["params"]["binning"]["metabat2"]["min_contig"],
        seed = config["params"]["binning"]["metabat2"]["seed"]
    shell:
        '''
        mkdir -p {params.bins_dir}
        metabat2 -i {input.asmfa} -a {input.depth} -o {params.bin_prefix} -m {params.min_contig} --seed {params.seed} -v > {log}
        touch {output.default}
        '''
