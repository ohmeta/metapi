def clean_reads(wildcards):
    if config["params"]["begin"] == "assembly":
        r1 = get_sample_id(_samples, wildcards, "fq1")
        r2 = get_sample_id(_samples, wildcards, "fq2")
        return [r1, r2]
    elif config["params"]["rmhost"]["do"]:
        return expand("{rmhost}/{sample}.rmhost.{read}.fq.gz",
                      rmhost=config["results"]["rmhost"],
                      sample=wildcards.sample,
                      read=["1", "2"])
    else:
        return expand("{trimming}/{sample}.trimmed.{read}.fq.gz",
                      trimming=config["results"]["trimming"],
                      sample=wildcards.sample,
                      read=["1", "2"])

rule metaquast:
    input:
        reads = clean_reads,
        scaftigs = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        report = os.path.join(config["results"]["metaquast"], "{sample}.{assembler}.metaquast_out/report.html"),
        icarus = os.path.join(config["results"]["metaquast"], "{sample}.{assembler}.metaquast_out/icarus.html"),
        combined_reference_tsv = os.path.join(config["results"]["metaquast"],
                                              "{sample}.{assembler}.metaquast_out/combined_reference/report.tsv"),
        icarus_viewers = directory(os.path.join(config["results"]["metaquast"],
                                                "{sample}.{assembler}.metaquast_out/icarus_viewers")),
        krona_charts = directory(os.path.join(config["results"]["metaquast"],
                                              "{sample}.{assembler}.metaquast_out/krona_charts")),
        not_aligned = directory(os.path.join(config["results"]["metaquast"],
                                             "{sample}.{assembler}.metaquast_out/not_aligned")),
        runs_per_reference = directory(os.path.join(config["results"]["metaquast"],
                                                    "{sample}.{assembler}.metaquast_out/runs_per_reference")),
        summary = directory(os.path.join(config["results"]["metaquast"],
                                         "{sample}.{assembler}.metaquast_out/summary"))
    log:
        os.path.join(config["logs"]["metaquast"], "{sample}.{assembler}.metaquast.log")
    params:
        output_dir = os.path.join(config["results"]["metaquast"], "{sample}.{assembler}.metaquast_out"),
        labels = "{sample}.{assembler}",
        metaquast_env = config["params"]["metaquast"]["env"]
    threads:
        config["params"]["metaquast"]["threads"]
    shell:
        '''
        set +u; source activate {params.metaquast_env}; set -u;
        metaquast.py {input.scaftigs} \
        --pe1 {input.reads[0]} \
        --pe2 {input.reads[1]} \
        -o {params.output_dir} \
        --labels {params.labels} \
        --threads {threads} 2> {log}
        '''

rule multiqc_metaquast:
    input:
        expand("{metaquast}/{sample}.{assembler}.metaquast_out/combined_reference/report.tsv",
               metaquast=config["results"]["metaquast"],
               assembler=config["params"]["assembler"],
               sample=_samples.index)
    output:
        html = os.path.join(config["results"]["metaquast"], "metaquast_multiqc_report.html"),
        data_dir = directory(os.path.join(config["results"]["metaquast"],
                                          "metaquast_multiqc_report_data"))
    log:
        os.path.join(config["logs"]["metaquast"], "multiqc_metaquast.log")
    params:
        outdir = config["results"]["metaquast"]
    shell:
        '''
        multiqc --outdir {params.outdir} --title metaquast --module quast {input} 2> {log}
        '''
