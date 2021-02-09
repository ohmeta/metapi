if config["params"]["binning"]["vamb"]["do"]:
    rule binning_vamb_combine_scaftigs:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{sample}.{{assembler}}.out/{sample}.{{assembler}}.scaftigs.fa.gz"),
                   sample=SAMPLES.index.unique())
        output:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.fa.gz")
        params:
            min_contig = config["params"]["binning"]["vamb"]["min_contig"]
        shell:
            '''
            concatenate.py {output} {input} -m {params.min_contig}
            '''


    rule binning_vamb_dict_scaftigs:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.fa.gz")
        output:
            os.path.join(
                config["output"]["alignment"],
                "index/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.dict")
        log:
            os.path.join(config["output"]["alignment"],
                         "logs/binning_vamb_dict_scaftigs_{assembler}.log")
        shell:
            '''
            samtools dict {input} | cut -f1-3 > {output} 2> {log}
            '''


    rule binning_vamb_index_scaftigs:
        input:
            os.path.join(
                config["output"]["assembly"],
                "scaftigs/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.fa.gz")
        output:
            os.path.join(
                config["output"]["alignment"],
                "index/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.minimap2.mmi")
        log:
            os.path.join(config["output"]["alignment"],
                         "logs/binning_vamb_index_scaftigs_{assembler}.log")
        params:
            index_size = config["params"]["binning"]["vamb"]["index_size"]
        shell:
            '''
            minimap2 -I {params.index_size} -d {output} {input} 2> {log}
            '''


    rule binning_vamb_align_scaftigs:
        input:
            reads = assembly_input_with_short_reads,
            index = os.path.join(
                config["output"]["alignment"],
                "index/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.minimap2.mmi"),
            dict = os.path.join(
                config["output"]["alignment"],
                "index/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.dict")
        output:
            flagstat = os.path.join(
                config["output"]["alignment"],
                "report/flagstat_minimap2/{sample}.{assembler}.align2combined_scaftigs.flagstat"),
            bam = os.path.join(
                config["output"]["alignment"],
                "bam/all.{assembler}.combined.out/{sample}.minimap2.out/{sample}.align2combined_scaftigs.sorted.bam") \
                if config["params"]["binning"]["vamb"]["save_bam"] else \
                   temp(os.path.join(
                       config["output"]["alignment"],
                       "bam/all.{assembler}.combined.out/{sample}.minimap2.out/{sample}.align2combined_scaftigs.sorted.bam"))
        log:
            os.path.join(config["output"]["alignment"],
                         "logs/alignment_multisplit/{sample}.{assembler}.align.reads2combined_scaftigs.log")
        threads:
            config["params"]["alignment"]["threads"]
        shell:
            '''
            rm -rf {output.bam}*

            minimap2 -t {threads} -ax sr {input.index} {input.reads} 2> {log} |
            tee >(samtools flagstat \
                  -@{threads} - > {output.flagstat}) | \
            grep -v "^@" | \
            cat {input.dict} - | \
            samtools view -F 3584 -b - |
            samtools sort -@{threads} -T {output.bam} -O BAM -o {output.bam} -
            '''


    rule binning_vamb_coverage:
        input:
            bam = os.path.join(
                config["output"]["alignment"],
                "bam/all.{assembler}.combined.out/{sample}.minimap2.out/{sample}.align2combined_scaftigs.sorted.bam")
        output:
            raw_jgi = os.path.join(
                config["output"]["multisplit_binning"],
                "coverage/{sample}.{assembler}.out/{sample}.align2combined_scaftigs.raw.jgi"),
            cut_jgi = os.path.join(
                config["output"]["multisplit_binning"],
                "coverage/{sample}.{assembler}.out/{sample}.align2combined_scaftigs.cut.jgi")
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/coverage/{sample}.{assembler}.align2combined_scaftigs.jgi.coverage.log")
        shell:
            '''
            jgi_summarize_bam_contig_depths \
            --noIntraDepthVariance --outputDepth {output.raw_jgi} {input.bam} 2> {log}

            cut -f1-3 --complement {output.raw_jgi} > {output.cut_jgi} 2>> {log}
            '''


    rule binning_vamb_gen_abundance_matrix:
        input:
            raw_jgi_first = os.path.join(
                config["output"]["multisplit_binning"],
                "coverage/%s.{assembler}.out/%s.align2combined_scaftigs.raw.jgi" % \
                (SAMPLES.index.unique()[0], SAMPLES.index.unique()[0])),
            cut_jgi = expand(os.path.join(
                config["output"]["multisplit_binning"],
                "coverage/{sample}.{{assembler}}.out/{sample}.align2combined_scaftigs.cut.jgi"),
                             sample=SAMPLES.index.unique())
        output:
            os.path.join(config["output"]["multisplit_binning"],
                         "matrix/all.{assembler}.align2combined_scaftigs.jgi.abundance.matrix.tsv")
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/binning_vamb_gen_abundance_matrix_{assembler}.log")
        shell:
            '''
            cut -f1-3 {input.raw_jgi_first} > {output}.column1to3
            paste {output}.column1to3 {input.cut_jgi} > {output} 2> {log}
            rm -rf {output}.column1to3
            '''


    rule binning_vamb:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/all.{assembler}.combined.out/all.{assembler}.combined.scaftigs.fa.gz"),
            matrix = os.path.join(
                config["output"]["multisplit_binning"],
                "matrix/all.{assembler}.align2combined_scaftigs.jgi.abundance.matrix.tsv")
        output:
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/clusters.tsv"),
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/latent.npz"),
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/lengths.npz"),
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/log.txt"),
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/model.pt"),
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/mask.npz"),
            os.path.join(config["output"]["multisplit_binning"],
                         "bins/all.{assembler}.combined.out/vamb/tnf.npz"),
            directory(os.path.join(config["output"]["multisplit_binning"],
                                   "bins/all.{assembler}.combined.out/vamb/bins"))
        log:
            os.path.join(config["output"]["multisplit_binning"],
                         "logs/binning/all.{assembler}.vamb.binning.log")
        params:
            outdir = os.path.join(
                config["output"]["multisplit_binning"],
                "bins/all.{assembler}.combined.out/vamb"),
            outdir_base = os.path.join(
                config["output"]["multisplit_binning"],
                "bins/all.{assembler}.combined.out/"),
            min_contig = config["params"]["binning"]["vamb"]["min_contig"]
        shell:
            '''
            rm -rf {params.outdir}
            mkdir -p {params.outdir_base}

            vamb \
            --outdir {params.outdir} \
            --fasta {input.scaftigs} \
            --jgi {input.matrix} \
            -o C -m {params.min_contig} --minfasta 500000 \
            2> {log}
            '''


    rule binning_vamb_all:
        input:
            expand(os.path.join(
                config["output"]["multisplit_binning"],
                "bins/all.{assembler}.combined.out/vamb/{results}"),
                   assembler=ASSEMBLERS,
                   results=["clusters.tsv", "latent.npz", "lengths.npz",
                            "log.txt", "model.pt", "mask.npz", "tnf.npz", "bins"])

else:
    rule binning_vamb_all:
        input:


rule multisplit_binning_all:
    input:
        rules.binning_vamb_all.input


rule binning_all:
    input:
        rules.single_binning_all.input,
        rules.cobinning_all.input,
        rules.multisplit_binning_all.input
