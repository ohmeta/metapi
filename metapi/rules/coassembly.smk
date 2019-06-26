def coassembly_inputs(read):
    if config["params"]["begin"] == "assembly":
        if config["params"]["reads_format"] == "fastq":
            if read == "1":
                return _samples.fq1
            elif read == "2":
                return _samples.fq2
        elif config["params"]["reads_format"] == "sra":
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
        os.path.join(config["logs"]["coassembly"]["megahit"], "megahit.coassembly.log")
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


rule demultiplex_kraken2_reads:
    input:
        reads = clean_reads,
        kraken2_output = os.path.join(config["results"]["classification"]["kraken2"], "{sample}.kraken2.output")
    output:
        done = os.path.join(config["results"]["coassembly"]["demultiplex_kraken2"], "{sample}.demultiplex_out/{sample}.demultiplex.done"),
        r1 = temp(os.path.join(config["results"]["coassembly"]["demultiplex_kraken2"], "{sample}.demultiplex_out/{sample}.1.fq")),
        r2 = temp(os.path.join(config["results"]["coassembly"]["demultiplex_kraken2"], "{sample}.demultiplex_out/{sample}.2.fq"))
    params:
        rank = config["params"]["coassembly"]["demultiplex_kraken2"]["rank"],
        taxadb = config["params"]["coassembly"]["demultiplex_kraken2"]["taxadb"],
        prefix = os.path.join(config["results"]["coassembly"]["demultiplex_kraken2"], "{sample}.demultiplex_out/{sample}"),
        change_seq_id = "--change_seq_id" if config["params"]["coassembly"]["demultiplex_kraken2"]["change_seq_id"] else ""
    log:
        os.path.join(config["logs"]["coassembly"]["demultiplex_kraken2"], "{sample}.demultiplex_kraken2.log")
    run:
        shell('''pigz -kdc {input.reads[0]} > {output.r1}''')
        shell('''pigz -kdc {input.reads[1]} > {output.r2}''')
        from metapi import demultiplexer
        demultiplexer.main(["--r1", output.r1,
                            "--r2", output.r2,
                            "--kraken2_output", input.kraken2_output,
                            "--rank", params.rank,
                            "--taxadb", params.taxadb,
                            "--prefix", params.prefix,
                            '--log', log,
                            params.change_seq_id])
        shell('''echo done > {output.done}''')


rule merge_kraken2_reads:
    input:
        expand("{demultiplex_kraken2}/{sample}.demultiplex_out/{sample}.demultiplex.done",
               demultiplex_kraken2=config["results"]["coassembly"]["demultiplex_kraken2"],
               sample=_samples.index.unique())
    output:
        os.path.join(config["results"]["coassembly"]["demultiplex_kraken2"], "merged")
    run:
        from glob import glob
        import os
        taxid_reads = {}
        for demultiplex_out in input:
            for r1 in glob(os.path.dirname(demultiplex_out) + "/*.1.fq.gz"):
                taxid = os.path.basename(r1).split(".")[-4]
                if taxid in taxid_reads:
                    taxid_reads[taxid]["r1"].append(r1)
                    taxid_reads[taxid]["r2"].append(r1.replace("1.fq", "2.fq"))
                else:
                    taxid_reads[taxid] = {}
                    taxid_reads[taxid]["r1"] = [r1]
                    taxid_reads[taxid]["r2"] = [r1.replace("1.fq", "2.fq")]
        
        for k,v in taxid_reads.items():
            r1_ = os.path.join(output, "%s.1.fq.gz" % k)
            r2_ = os.path.join(output, "%s.2.fq.gz" % k)
            r1_str = " ".join(v["r1"])
            r2_str = " ".join(v["r2"])
            shell('''cat %s > %s''' % (r1_str, r1_))
            shell('''cat %s > %s''' % (r2_str, r2_))
            shell('''rm -rf %s''' % r1_str)
            shell('''rm -rf %s''' % r2_str)
