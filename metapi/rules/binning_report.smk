if len(BINNERS_CHECKM) != 0:
    rule binning_report:
        input:
            lambda wildcards: get_binning_done(wildcards, [wildcards.binner_checkm])
        output:
            directory(
                os.path.join(
                    config["output"]["binning"],
                    "report/{assembler}_{binner_checkm}_stats/{binning_group}.{assembly_group}"))
        params:
            binning_group = "{binning_group}",
            assembly_group = "{assembly_group}",
            assembler = "{assembler}",
            binner = "{binner_checkm}"
        priority:
            35
        run:
            import glob

            shell('''rm -rf {output}''')
            shell('''mkdir -p {output}''')

            bin_list =  glob.glob(os.path.dirname(input[0]) + "/*.fa.gz")
            header_list = [
                "binning_group", "assembly_group", "bin_id", "bin_file", "assembler", "binner",
                "chr", "length", "#A", "#C", "#G", "#T",
                "#2", "#3", "#4", "#CpG", "#tv", "#ts", "#CpG-ts"]
            header_name = "\\t".join(header_list)

            for bin_fa in bin_list:
                bin_id = os.path.basename(os.path.splitext(os.path.splitext(bin_fa)[0])[0])
                bin_file = os.path.abspath(bin_fa)
                header_content = "\\t".join([params.binning_group, params.assembly_group, bin_id, bin_file, params.assembler, params.binner])
                stats_file = os.path.join(output[0], bin_id + ".seqtk.comp.tsv.gz")

                shell(
                    '''
                    seqtk comp %s | \
                    awk \
                    'BEGIN \
                    {{print "%s"}}; \
                    {{print "%s" "\t" $0}}' | \
                    gzip -c > %s
                    ''' % (bin_fa, header_name, header_content, stats_file))


    rule binning_report_merge:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "report/{{assembler}}_{{binner_checkm}}_stats/{binning_group}.{assembly_group}"),
                zip,
                binning_group=ASSEMBLY_GROUP["binning_group"],
                assembly_group=ASSEMBLY_GROUP["assembly_group"])
        output:
            summary = os.path.join(
                config["output"]["binning"],
                "report/assembly_stats_{assembler}_{binner_checkm}.tsv.gz")
        params:
            min_length = config["params"]["assembly"]["report"]["min_length"],
            len_ranges = config["params"]["assembly"]["report"]["len_ranges"]
        threads:
            config["params"]["binning"]["threads"]
        run:
            import glob
            comp_list = []
            for i in input:
                comp_list += glob.glob(i + "/*.seqtk.comp.tsv.gz")

            if len(comp_list) != 0:
                metapi.assembler_init(
                    params.len_ranges,
                    ["binning_group", "assembly_group", "bin_id", "bin_file", "assembler", "binner"])
                comp_list_ = [(j, params.min_length) for j in comp_list]
                metapi.merge(
                    comp_list_, metapi.parse_assembly,
                    threads, output=output.summary)
            else:
                shell('''touch {output.summary}''')


    rule binning_report_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "report/assembly_stats_{assembler}_{binner_checkm}.tsv.gz"),
                assembler=ASSEMBLERS,
                binner_checkm=BINNERS_CHECKM)

else:
    rule binning_report_all:
        input:


rule binning_all:
    input:
        rules.binning_metabat2_all.input,
        rules.binning_maxbin2_all.input,
        rules.binning_concoct_all.input,
        rules.binning_graphbin2_all.input,
        rules.binning_vamb_all.input,
        rules.binning_dastools_all.input,
        rules.binning_report_all.input


localrules:
    binning_report_all,
    binning_all