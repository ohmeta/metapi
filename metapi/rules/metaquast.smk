rule metaquast_megahit:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    output:
        report = os.path.join(config["results"]["metaquast"], "{sample}.metaquast_out/report.html"),
        icarus = os.path.join(config["results"]["metaquast"], "{sample}.metaquast_out/icarus.html"),
        combined_reference = directory(os.path.join(config["results"]["metaquast"],
                                                    "{sample}.metaquast_out/combined_reference")),
        icarus_viewers = directory(os.path.join(config["results"]["metaquast"],
                                                "{sample}.metaquast_out/icarus_viewers")),
        krona_charts = directory(os.path.join(config["results"]["metaquast"],
                                              "{sample}.metaquast_out/krona_charts")),
        not_aligned = directory(os.path.join(config["results"]["metaquast"],
                                             "{sample}.metaquast_out/not_aligned")),
        runs_per_reference = directory(os.path.join(config["results"]["metaquast"],
                                                    "{sample}.metaquast_out/runs_per_reference")),
        summary = directory(os.path.join(config["results"]["metaquast"],
                                         "{sample}.metaquast_out/summary"))
    log:
        os.path.join(config["logs"]["metaquast"], "{sample}.metaquast.log")
    params:
        output_dir = os.path.join(config["results"]["metaquast"], "{sample}.metaquast_out"),
        min_contig = config["params"]["metaquast"]["min_contig"]
    threads:
        config["params"]["metaquast"]["threads"]
    shell:
        '''
        metaquast.py {input} -o {params.output_dir} \
        --min-contig {params.min_contig} \
        --threads {threads} 2> {log}
        '''

'''
rule metaquast_idba_ud:
    input:
    output:
    log:
    params:
    threads:
    shell:

rule metaquast_metaspades:
    input:
    output:
    log:
    params:
    threads:
    shell:
'''
