def clean_reads(wildcards):
    if config["params"]["begin"] == "assembly":
        if config["params"]["type"] == "fastq":
            r1 = get_sample_id(_samples, wildcards, "fq1")
            r2 = get_sample_id(_samples, wildcards, "fq2")
            return [r1, r2]
        elif config["params"]["type"] == "sra":
            return expand("{sra2fq}/{sample}.{read}.fq.gz",
                          sra2fq=config["results"]["sra2fq"],
                          sample=wildcards.sample,
                          read=[1, 2])
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

rule assembly_megahit:
    input:
        reads = clean_reads
    output:
        contigs = protected(os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.megahit.scaftigs.fa.gz")),
        temp_file = temp(directory(os.path.join(config["results"]["assembly"], "{sample}.megahit_out/intermediate_contigs")))
    params:
        min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.megahit_out"),
        contigs = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.{megahit}.contigs.fa"),
        out_prefix = "{sample}.megahit"
    threads:
        config["params"]["assembly"]["megahit"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["megahit"], "{sample}.megahit.log")
    shell:
        '''
        rm -rf {params.out_dir}
        megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} --min-contig-len {params.min_contig} --out-dir {params.out_dir} --out-prefix {params.out_prefix} 2> {log}
        sed -i 's#^>#>'"{params.out_prefix}"'_#g' {params.contigs}
        pigz {params.contigs}
        mv {params.contigs}.gz {output.contigs}
        '''

rule assembly_idba_ud:
    input:
        reads = clean_reads
    output:
        protected(os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out/{sample}.idba_ud.scaftigs.fa.gz"))
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
        if [[ $(file --mime-type -b {input.reads[0]}) == application/x-gzip ]]; then
            zcat {input.reads[0]} > {params.r1}
            zcat {input.reads[1]} > {params.r2}
        else
            ln -s $(realpath {input.reads[0]}) {params.r1}
            ln -s $(realpath {input.reads[0]}) {params.r2}
        fi
        fq2fa --merge {params.r1} {params.r2} {params.pe_fa}
        idba_ud -r {params.pe_fa} --mink {params.mink} --maxk {params.maxk} --step {params.step} --min_contig {params.min_contig} -o {params.out_dir} --num_threads {threads} --pre_correction 2> {log}
        pigz {params.out_dir}/scaffold.fa
        mv {params.out_dir}/scaffold.fa.gz {output}
        rm -rf {params.out_dir}/kmer {params.out_dir}/contig-* {params.out_dir}/align-* {params.out_dir}/graph-* {params.out_dir}/local-contig-*
        '''

kmer_list = ["21", "33", "55"]
if len(config["params"]["assembly"]["metaspades"]["kmers"]) > 0:
    kmer_list = config["params"]["assembly"]["metaspades"]["kmers"]

def get_kmer_dirs(wildcards):
    return expand(directory(temp(os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/K{kmer}"))),
                  sample=wildcards.sample, kmer=kmer_list)

rule assembly_metaspades:
    input:
        reads = clean_reads
    output:
        scaftigs = protected(os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/{sample}.metaspades.scaftigs.fa.gz"))
    params:
        kmers = "auto" if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out"),
        corrected = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/corrected"),
        kmer_dirs = get_kmer_dirs,
        only_assembler = "--only-assembler" if config["params"]["assembly"]["metaspades"]["only_assembler"] else "",
        only_save_scaftigs = "true" if config["params"]["assembly"]["metaspades"]["only_save_scaftigs"] else "false",
        tar_results = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/{sample}.metaspades.tar")
    threads:
        config["params"]["assembly"]["metaspades"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["metaspades"], "{sample}.metaspades.log")
    shell:
        '''
        metaspades.py \
        -1 {input.reads[0]} \
        -2 {input.reads[1]} \
        -k {params.kmers} \
        {params.only_assembler} \
        --threads {threads} \
        -o {params.out_dir} 2> {log}
        pigz -p {threads} {params.out_dir}/scaffolds.fasta
        mv {params.out_dir}/scaffolds.fasta.gz {output.scaftigs}
        rm -rf {params.kmer_dirs}
        rm -rf {params.corrected}
        
        if {params.only_save_scaftigs}; then
            find {params.out_dir} -type f ! -wholename "{output.scaftigs}" -delete
        else
            find {params.out_dir} -type f ! -wholename "{output.scaftigs}" ! -wholename "{params.tar_results}" | xargs -I % sh -c 'tar -rf {params.tar_results} %; rm -rf %'
            pigz -p {threads} {params.tar_results}
        fi
        
        rm -rf {params.out_dir}/tmp
        rm -rf {params.out_dir}/misc
        '''
