rule quast_megahit:
    input:
        os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    output:
        report = os.path.join(config["results"]["quast"], "{sample}.metaquast_out/report.html"),
        icarus = os.path.join(config["results"]["quast"],
                              "{sample}.quast_out/icarus.html"),
        contigs_reports = directory(os.path.join(config["results"]["quast"],
                                                 "{sample}.metaquast_out/contigs_reports")),
        k_mer_stats = directory(os.path.join(config["results"]["quast"],
                                             "{sample}.metaquast_out/k_mer_stats")),
        reads_stats = directory(os.path.join(config["results"]["quast"],
                                             "{sample}.metaquast_out/reads_stats"))
    log:
        config["logs"]["quast"]
    params:
        output_dir = os.path.join(config["results"]["quast"], "{sample}.quast_out"),
        min_contig = config["params"]["quast"]["min_contig"]
    threads:
        config["params"]["quast"]["threads"]
    shell:
        '''
        metaquast.py {input} -o {params.output_dir} \
        --min-contig {params.min_contig} \
        --threads {threads} \
        '''

'''
rule quast_idba_ud:
    input:
    output:
    log:
    params:
    threads:
    shell:

rule quast_metaspades:
    input:
    output:
    log:
    params:
    threads:
    shell:
'''
