def assembly_inputs(wildcards):
    if config["params"]["rmhost"]["do"]:
        return expand("{rmhost}/{sample}.rmhost.{read}.fq.gz",
                      rmhost=config["results"]["rmhost"],
                      sample=wildcards.sample,
                      read=["1", "2"])
    else:
        return expand("{trimming}/{sample}.trimmed.{read}.fq.gz",
                      trimming=config["results"]["trimming"],
                      sample=wildcards.sample,
                      read=["1", "2"])

rule assembly_megahit:
    input:
        reads = assembly_inputs
    output:
        os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.contigs.fa.gz")
    params:
        min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.megahit_out"),
        out_prefix = "{sample}"
    threads:
        config["params"]["assembly"]["megahit"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["megahit"], "{sample}.megahit.log")
    shell:
        '''
        rm -rf {params.out_dir}
        megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} --min-contig-len {params.min_contig} --out-dir {params.out_dir} --out-prefix {params.out_prefix} 2> {log}
        pigz {params.out_dir}/{params.out_prefix}.contigs.fa
        '''

rule assembly_idba_ud:
    input:
        reads = assembly_inputs
    output:
        os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out/{sample}.scaffolds.fa.gz")
    params:
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out"),
        r1 = temp(os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out/{sample}.r1.fq")),
        r2 = temp(os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out/{sample}.r2.fq")),
        pe_fa = temp(os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out/{sample}.pe.fa")),
        mink = config["params"]["assembly"]["idba_ud"]["mink"],
        maxk = config["params"]["assembly"]["idba_ud"]["maxk"],
        step = config["params"]["assembly"]["idba_ud"]["step"],
        min_contig = config["params"]["assembly"]["idba_ud"]["min_contig"]
    threads:
        config["params"]["assembly"]["idba_ud"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["idba_ud"], "{sample}.idba_ud.log")
    shell:
        '''
        rm -rf {params.out_dir}
        mkdir {params.out_dir}
        zcat {input.reads[0]} > {params.r1}
        zcat {input.reads[1]} > {params.r2}
        fq2fa --merge {params.r1} {params.r2} {params.pe_fa}
        idba_ud -r {params.pe_fa} --mink {params.mink} --maxk {params.maxk} --step {params.step} --min_contig {params.min_contig} -o {params.out_dir} --num_threads {threads} --pre_correction 2> {log}
        pigz {params.out_dir}/scaffold.fa
        mv {params.out_dir}/scaffold.fa.gz {output}
        rm -rf {params.out_dir}/kmer {params.out_dir}/contig-* {params.out_dir}/align-* {params.out_dir}/graph-* {params.out_dir}/local-contig-*
        '''

rule assembly_metaspades:
    input:
        reads = assembly_inputs
    output:
        os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/{sample}.scaffolds.fa.gz")
    params:
        memory = config["params"]["assembly"]["metaspades"]["memory"],
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out")
    threads:
        config["params"]["assembly"]["metaspades"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["metaspades"], "{sample}.metaspades.log")
    shell:
        '''
        metaspades.py -1 {input.reads[0]} -2 {input.reads[1]} --threads {threads} --memory {params.memory} -o {params.out_dir} 2> {log}
        pigz {params.out_dir}/scaffolds.fasta
        mv {params.out_dir}/scaffolds.fasta.gz {output}
        '''
