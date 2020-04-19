if config["params"]["binning"]["metabat2"]["do"]:
    rule binning_metabat2_coverage:
        input:
            bam = os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam"),
            bai = os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam.bai")
        output:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.metabat2.coverage")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{sample}.{assembler}.metabat2.coverage.log")
        params:
            output_dir = os.path.join(config["output"]["binning"],
                                      "coverage/{sample}.{assembler}.out")
        shell:
            '''
            mkdir -p {params.output_dir}
            jgi_summarize_bam_contig_depths \
            --outputDepth {output.coverage} \
            {input.bam} \
            2> {log}
            '''


    rule binning_metabat2:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz"),
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.metabat2.coverage")
        output:
            bins_dir = directory(os.path.join(config["output"]["binning"],
                                              "bins/{sample}.{assembler}.out/metabat2"))
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{sample}.{assembler}.metabat2.binning.log")
        params:
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "bins/{sample}.{assembler}.out/metabat2/{sample}.{assembler}.metabat2.bin"),
            min_contig = config["params"]["binning"]["metabat2"]["min_contig"],
            seed = config["params"]["binning"]["metabat2"]["seed"]
        shell:
            '''
            rm -rf {output.bins_dir}
            mkdir -p {output.bins_dir}

            metabat2 \
            -i {input.scaftigs} \
            -a {input.coverage} \
            -o {params.bin_prefix} \
            -m {params.min_contig} \
            --seed {params.seed} -v \
            > {log}
            '''


    rule binning_metabat2_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins/{sample}.{assembler}.out/metabat2"),
                assembler=ASSEMBLERS,
                sample=SAMPLES.index.unique())


if config["params"]["binning"]["maxbin2"]["do"]:
    rule binning_maxbin2_coverage:
        input:
            bam = os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam"),
            bai = os.path.join(
                config["output"]["alignment"],
                "bam/{sample}.{assembler}.out/{sample}.{assembler}.align2scaftigs.sorted.bam.bai")
        output:
            coverage_bb = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.bbmap.coverage"),
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.maxbin2.coverage")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{sample}.{assembler}.maxbin2.coverage.log")
        params:
            output_dir = os.path.join(config["output"]["binning"],
                                      "coverage/{sample}.{assembler}.out")
        shell:
            '''
            mkdir -p {params.output_dir}
            pileup.sh in={input.bam} out={output.coverage_bb} 2> {log}
            awk '{print $1 "\t" $5}' {output.coverage_bb} | grep -v '^#' > {output.coverage}
            '''


    rule binning_maxbin2:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{assembler}.out/{sample}.{assembler}.scaftigs.fa.gz"),
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.maxbin2.coverage")
        output:
            bins_dir = directory(os.path.join(config["output"]["binning"],
                                              "bins/{sample}.{assembler}.out/maxbin2"))
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{sample}.{assembler}.maxbin2.binning.log")
        params:
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "bins/{sample}.{assembler}.out/maxbin2/{sample}.{assembler}.maxbin2.bin")
        threads:
            config["params"]["binning"]["maxbin2"]["threads"]
        shell:
            '''
            rm -rf {output.bins_dir}
            mkdir -p {output.bins_dir}

            run_MaxBin.pl \
            -thread {threads} \
            -contig {input.scaftigs} \
            -abund {input.coverage} \
            -out {params.bin_prefix} \
            2> {log}
            '''


    rule binning_maxbin2_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins/{sample}.{assembler}.out/maxbin2"),
                assembler=ASSEMBLERS,
                sample=SAMPLES.index.unique())


if len(BINNERS) != 0:
    rule binning_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins/{sample}.{assembler}.out/{binner}"),
                assembler=ASSEMBLERS,
                binner=BINNERS,
                sample=SAMPLES.index.unique())
else:
    rule binning_all:
        input:
