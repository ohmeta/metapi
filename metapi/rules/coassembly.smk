def coassembly_inputs(wildcards):
    if config["params"]["rmhost"]["do"]:
        r1 = expand("{rmhost}/{sample}.rmhost.1.fq.gz",
                       rmhost=config["results"]["rmhost"],
                       sample=_samples.index)
        r2 = expand("{rmhost}/{sample}.rmhost.2.fq.gz",
                       rmhost=config["results"]["rmhost"],
                       sample=_samples.index)
        r1_ = ",".join(r1)
        r2_ = ",".join(r2)
        return [r1_, r2_]
    else:
        r1 = expand("{trimming}/{sample}.trimmed.1.fq.gz",
                       trimming=config["results"]["trimming"],
                       sample=_samples.index)
        r2 = expand("{trimming}/{sample}.trimmed.2.fq.gz",
                       trimming=config["results"]["trimming"],
                       sample=_samples.index)
        r1_ = ",".join(r1)
        r2_ = ",".join(r2)
        return [r1_, r2_]


rule coassembly_megahit:
    input:
        coassembly_inputs
    output:
        contigs = os.path.join(config["results"]["coassembly"]["megahit"], "contigs.fa.gz"),
        temp_file = temp(directory(os.path.join(config["results"]["coassembly"]["megahit"],
                                                "intermediate_contigs")))
    log:
        os.path.join(config["logs"]["coassembly"], "megahit.coassembly.log")
    params:
        min_contig = config["params"]["coassembly"]["megahit"]["min_contig"],
        out_dir = config["results"]["coassembly"]["megahit"]
    threads:
        config["params"]["coassembly"]["megahit"]["threads"]
    shell:
        '''
        meghait -1 {input[0]} -2 {input[1]} -t {threads} --min-contig-len {params.min_contig} --out-dir {params.out_dir} 2> {log}
        '''
