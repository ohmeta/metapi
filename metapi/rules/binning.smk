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
        priority:
            30
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{sample}.{assembler}.metabat2.coverage.log")
        params:
            output_dir = os.path.join(config["output"]["binning"],
                                      "coverage/{sample}.{assembler}.out")
        shell:
            '''
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
        priority:
            30
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

else:
    rule binning_metabat2_all:
        input:


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
        priority:
            30
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{sample}.{assembler}.maxbin2.coverage.log")
        params:
            output_dir = os.path.join(config["output"]["binning"],
                                      "coverage/{sample}.{assembler}.out")
        shell:
            '''
            pileup.sh in={input.bam} out={output.coverage_bb} 2> {log}
            awk '{{print $1 "\t" $5}}' {output.coverage_bb} | grep -v '^#' > {output.coverage}
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
        priority:
            30
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{sample}.{assembler}.maxbin2.binning.log")
        params:
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "bins/{sample}.{assembler}.out/maxbin2/{sample}.{assembler}.maxbin2.bin")
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
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

else:
    rule binning_maxbin2_all:
        input:


if len(BINNERS) != 0:
    rule binning_report:
        input:
            bins_dir = os.path.join(
                config["output"]["binning"],
                "bins/{sample}.{assembler}.out/{binner}")
        output:
            report_dir = directory(
                os.path.join(
                    config["output"]["binning"],
                    "report/{assembler}_{binner}_stats/{sample}"))
        priority:
            35
        params:
            sample_id = "{sample}",
            assembler = "{assembler}",
            binner = "{binner}"
        run:
            import glob

            shell('''rm -rf {output.report_dir}''')
            shell('''mkdir -p {output.report_dir}''')

            bin_list =  glob.glob(input.bins_dir + "/*bin*fa")
            header_list = ["sample_id", "bin_id", "assembler", "binner",
                           "chr", "length", "#A", "#C", "#G", "#T",
                           "#2", "#3", "#4", "#CpG", "#tv", "#ts", "#CpG-ts"]
            header = "\\t".join(header_list)

            for bin_fa in bin_list:
                bin_id = os.path.basename(os.path.splitext(bin_fa)[0])
                header_ = "\\t".join([params.sample_id, bin_id,
                                      params.assembler, params.binner])
                stats_file = os.path.join(output.report_dir,
                                          bin_id + ".seqtk.comp.tsv.gz")

                shell(
                    '''
                    seqtk comp %s | \
                    awk \
                    'BEGIN \
                    {{print "%s"}}; \
                    {{print "%s" "\t" $0}}' | \
                    gzip -c > %s
                    ''' % (bin_fa, header, header_, stats_file))


    rule binning_report_merge:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "report/{{assembler}}_{{binner}}_stats/{sample}"),
                   sample=SAMPLES.index.unique())
        output:
            summary = os.path.join(
                config["output"]["binning"],
                "report/assembly_stats_{assembler}_{binner}.tsv")
        params:
            len_ranges = config["params"]["assembly"]["report"]["len_ranges"]
        threads:
            config["params"]["binning"]["threads"]
        run:
            import glob
            comp_list = []
            for i in input:
                comp_list += glob.glob(i + "/*bin*.seqtk.comp.tsv.gz")

            if len(comp_list) != 0:
                metapi.assembler_init(params.len_ranges,
                                      ["sample_id", "bin_id", "assembler", "binner"])
                metapi.merge(comp_list, metapi.parse_assembly,
                             threads, save=True,  output=output.summary)
            else:
                shell('''touch {output.summary}''')


    rule binning_report_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "report/assembly_stats_{assembler}_{binner}.tsv"),
                   assembler=ASSEMBLERS,
                   binner=BINNERS)

else:
    rule binning_report_all:
        input:


rule binning_all:
    input:
        rules.binning_metabat2_all.input,
        rules.binning_maxbin2_all.input,
        rules.binning_report_all.input,

        rules.alignment_all.input
