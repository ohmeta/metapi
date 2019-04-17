def coassembly_inputs(read):
    if config["params"]["begin"] == "assembly":
        if config["params"]["type"] == "fastq":
            if read == "1":
                return _samples.fq1
            elif read == "2":
                return _samples.fq2
        elif config["params"]["type"] == "sra":
            return expand("{sra2fq}/{sample}.{read}.fq.gz",
                          sra2fq=config["results"]["sra2fq"],
                          read=read,
                          sample=_samples.index.unique())
    elif config["params"]["rmhost"]["do"]:
        return expand("{rmhost}/{sample}.rmhost.{read}.fq.gz",
                      rmhost=config["results"]["rmhost"],
                      read=read,
                      sample=_samples.index.unique())
    else:
        return expand("{trimming}/{sample}.trimmed.{read}.fq.gz",
                      trimming=config["results"]["trimming"],
                      read=read,
                      sample=_samples.index.unique())

rule coassembly_megahit:
    input:
        r1 = coassembly_inputs("1"),
        r2 = coassembly_inputs("2")
    output:
        contigs = os.path.join(config["results"]["coassembly"]["megahit"], "final.contigs.fa.gz"),
        temp_file = temp(directory(os.path.join(config["results"]["coassembly"]["megahit"],
                                                "intermediate_contigs")))
    log:
        os.path.join(config["logs"]["coassembly"], "megahit.coassembly.log")
    params:
        min_contig = config["params"]["coassembly"]["megahit"]["min_contig"],
        out_dir = config["results"]["coassembly"]["megahit"]
    threads:
        config["params"]["coassembly"]["megahit"]["threads"]
    run:
        r1_str = ",".join(input.r1)
        r2_str = ",".join(input.r2)
        shell('''rm -rf {params.out_dir}
        megahit -1 {r1_str} -2 {r2_str} -t {threads} --min-contig-len {params.min_contig} --out-dir {params.out_dir} 2> {log}
        pigz {params.out_dir}/final.contigs.fa
        ''')
