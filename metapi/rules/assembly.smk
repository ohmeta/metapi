def clean_reads(wildcards):
    if config["params"]["begin"] == "assembly":
        if config["params"]["reads_format"] == "fastq":
            if IS_PE:
                return [sample.get_sample_id(_samples, wildcards, "fq1"),
                        sample.get_sample_id(_samples, wildcards, "fq2")]
            else:
                return [sample.get_sample_id(_samples, wildcards, "fq1")]
        elif config["params"]["reads_format"] == "sra":
            return expand("{sra2fq}/{sample}{read}.fq.gz",
                          sra2fq=config["results"]["sra2fq"],
                          sample=wildcards.sample,
                          read=[".1", ".2"] if IS_PE else "")
    elif config["params"]["rmhost"]["do"]:
        return expand("{rmhost}/{sample}.rmhost{read}.fq.gz",
                      rmhost=config["results"]["rmhost"],
                      sample=wildcards.sample,
                      read=[".1", ".2"] if IS_PE else "")
    else:
        return expand("{trimming}/{sample}.trimmed{read}.fq.gz",
                      trimming=config["results"]["trimming"],
                      sample=wildcards.sample,
                      read=[".1", ".2"] if IS_PE else "")

rule assembly_megahit:
    input:
        reads = clean_reads
    output:
        contigs = protected(os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.megahit.scaftigs.fa.gz")),
        temp_file = temp(directory(os.path.join(config["results"]["assembly"], "{sample}.megahit_out/intermediate_contigs")))
    params:
        min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.megahit_out"),
        contigs = os.path.join(config["results"]["assembly"], "{sample}.megahit_out/{sample}.megahit.contigs.fa"),
        out_prefix = "{sample}.megahit"
    threads:
        config["params"]["assembly"]["megahit"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["megahit"], "{sample}.megahit.log")
    run:
        shell("rm -rf {params.out_dir}")
        if IS_PE:
            shell("megahit -1 {input.reads[0]} -2 {input.reads[1]} -t {threads} --min-contig-len {params.min_contig} --out-dir {params.out_dir} --out-prefix {params.out_prefix} 2> {log}")
        else:
            shell("megahit -r {input.reads[0]} -t {threads} --min-contig-len {params.min_contig} --out-dir {params.out_dir} --out-prefix {params.out_prefix} 2> {log}")
        shell('''sed -i 's#^>#>'"{params.out_prefix}"'_#g' {params.contigs}''')
        shell("pigz {params.contigs}")
        shell("mv {params.contigs}.gz {output.contigs}")

rule assembly_idba_ud:
    input:
        reads = clean_reads
    output:
        protected(os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out/{sample}.idba_ud.scaftigs.fa.gz"))
    params:
        prefix = "{sample}",
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.idba_ud_out"),
        mink = config["params"]["assembly"]["idba_ud"]["mink"],
        maxk = config["params"]["assembly"]["idba_ud"]["maxk"],
        step = config["params"]["assembly"]["idba_ud"]["step"],
        min_contig = config["params"]["assembly"]["idba_ud"]["min_contig"]
    threads:
        config["params"]["assembly"]["idba_ud"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["idba_ud"], "{sample}.idba_ud.log")
    run:
        shell("rm -rf {params.out_dir}")
        shell("mkdir {params.out_dir}")
        reads = os.path.join(config["results"]["assembly"], "%s.idba_ud_out/%s.fa" % (params.prefix, params.prefix))
        if IS_PE:
            shell("seqtk mergepe {input.reads[0]} {input.reads[1]} | seqtk seq -A - > %s" % reads)
        else:
            shell("seqtk seq -A {intput.reads[0]} > %s" % reads)

        shell("idba_ud -r %s --mink {params.mink} --maxk {params.maxk} --step {params.step} \
              --min_contig {params.min_contig} \
              -o {params.out_dir} --num_threads {threads} --pre_correction 2> {log}" % reads)

        shell("pigz {params.out_dir}/scaffold.fa")
        shell("mv {params.out_dir}/scaffold.fa.gz {output}")
        shell("rm -rf %s" % reads)
        shell("rm -rf {params.out_dir}/kmer")
        shell("rm -rf {params.out_dir}/contig-*")
        shell("rm -rf {params.out_dir}/align-*")
        shell("rm -rf {params.out_dir}/graph-*")
        shell("rm -rf {params.out_dir}/local-contig-*")


metaspades_kmer_list = ["21", "33", "55"]
if len(config["params"]["assembly"]["metaspades"]["kmers"]) > 0:
    metaspades_kmer_list = config["params"]["assembly"]["metaspades"]["kmers"]

spades_kmer_list = ["21", "33", "55"]
if len(config["params"]["assembly"]["spades"]["kmers"]) > 0:
    spades_kmer_list = config["params"]["assembly"]["spades"]["kmers"]

def get_kmer_dirs(wildcards, asmer, kmer_list):
    return directory(temp(expand(os.path.join(config["results"]["assembly"],
                                              "{sample}.{asmer}_out/K{kmer}"),
                  asmer=asmer,
                  sample=wildcards.sample,
                  kmer=kmer_list)))

rule assembly_metaspades:
    input:
        reads = clean_reads
    output:
        scaftigs = protected(os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/{sample}.metaspades.scaftigs.fa.gz"))
    params:
        kmers = "auto" if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out"),
        corrected = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/corrected"),
        kmer_dirs = lambda wildcards: get_kmer_dirs(wildcards, "metaspades", metaspades_kmer_list),
        only_assembler = "--only-assembler" if config["params"]["assembly"]["metaspades"]["only_assembler"] else "",
        only_save_scaftigs = "true" if config["params"]["assembly"]["metaspades"]["only_save_scaftigs"] else "false",
        tar_results = os.path.join(config["results"]["assembly"], "{sample}.metaspades_out/{sample}.metaspades.tar")
    threads:
        config["params"]["assembly"]["metaspades"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["metaspades"], "{sample}.metaspades.log")
    run:
        if IS_PE:
            shell('''metaspades.py \
                  -1 {input.reads[0]} \
                  -2 {input.reads[1]} \
                  -k {params.kmers} \
                  {params.only_assembler} \
                  --threads {threads} \
                  -o {params.out_dir} 2> {log}''')
            shell("pigz -p {threads} {params.out_dir}/scaffolds.fasta")
            shell("mv {params.out_dir}/scaffolds.fasta.gz {output.scaftigs}")
            shell("rm -rf {params.kmer_dirs}")
            shell("rm -rf {params.corrected}")

            if config["params"]["assembly"]["metaspades"]["only_save_scaftigs"]:
                shell('''find {params.out_dir} -type f ! -wholename "{output.scaftigs}" -delete''')
            else:
                shell("""find {params.out_dir} -type f \
                      ! -wholename "{output.scaftigs}" \
                      ! -wholename "{params.tar_results}" \
                      | xargs -I % sh -c 'tar -rf {params.tar_results} %; rm -rf %'""")
                shell("pigz -p {threads} {params.tar_results}")

            shell("rm -rf {params.out_dir}/tmp")
            shell("rm -rf {params.out_dir}/misc")
        else:
            print("don't support single-end reads assembly using metaspades, you can try spades or megahit, idba_ud")
            sys.exit(1)

rule assembly_spades:
    input:
        reads = clean_reads
    output:
        scaftigs = protected(os.path.join(config["results"]["assembly"], "{sample}.spades_out/{sample}.spades.scaftigs.fa.gz"))
    params:
        kmers = "auto" if len(config["params"]["assembly"]["spades"]["kmers"]) == 0 else ",".join(config["params"]["assembly"]["spades"]["kmers"]),
        out_dir = os.path.join(config["results"]["assembly"], "{sample}.spades_out"),
        corrected = os.path.join(config["results"]["assembly"], "{sample}.spades_out/corrected"),
        kmer_dirs = lambda wildcards: get_kmer_dirs(wildcards, "spades", spades_kmer_list),
        only_assembler = "--only-assembler" if config["params"]["assembly"]["spades"]["only_assembler"] else "",
        only_save_scaftigs = "true" if config["params"]["assembly"]["spades"]["only_save_scaftigs"] else "false",
        tar_results = os.path.join(config["results"]["assembly"], "{sample}.spades_out/{sample}.spades.tar")
    threads:
        config["params"]["assembly"]["spades"]["threads"]
    log:
        os.path.join(config["logs"]["assembly"]["spades"], "{sample}.spades.log")
    run:
        if IS_PE:
            shell('''spades.py \
                  -1 {input.reads[0]} \
                  -2 {input.reads[1]} \
                  -k {params.kmers} \
                  {params.only_assembler} \
                  --threads {threads} \
                  -o {params.out_dir} 2> {log}''')
        else:
            shell('''spades.py \
                  -s {input.reads[0]} \
                  -k {params.kmers} \
                  {params.only_assembler} \
                  --threads {threads} \
                  -o {params.out_dir} 2> {log}''')

        shell("pigz -p {threads} {params.out_dir}/scaffolds.fasta")
        shell("mv {params.out_dir}/scaffolds.fasta.gz {output.scaftigs}")
        shell("rm -rf {params.kmer_dirs}")
        shell("rm -rf {params.corrected}")

        if config["params"]["assembly"]["spades"]["only_save_scaftigs"]:
            shell('''find {params.out_dir} -type f ! -wholename "{output.scaftigs}" -delete''')
        else:
            shell("""find {params.out_dir} -type f \
                  ! -wholename "{output.scaftigs}" \
                  ! -wholename "{params.tar_results}" \
                  | xargs -I % sh -c 'tar -rf {params.tar_results} %; rm -rf %'""")
            shell("pigz -p {threads} {params.tar_results}")

        shell("rm -rf {params.out_dir}/tmp")
        shell("rm -rf {params.out_dir}/misc")


rule assembly_report:
    input:
        scaftigs = os.path.join(config["results"]["assembly"],
                                "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz")
    output:
        report = os.path.join(config["results"]["assembly"],
                              "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.seqtk.comp.tsv.gz")
    params:
        sample_id = "{sample}"
    shell:
        """
        seqtk comp {input.scaftigs} | \
        awk 'BEGIN {{print "sample_id\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts"}}; \
             {{print "{params.sample_id}" "\t" $0}}' | \
        gzip -c > {output.report}
        """
