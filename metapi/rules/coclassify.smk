if config["params"]["classify"]["gtdbtk"]["do"]:
    rule coclassify_hmq_bins_gtdbtk:
        input:
            os.path.join(
                config["output"]["cocheckm"],
                "bins_hmq/{assembler_co}.{binner_checkm}.links")
        output:
            os.path.join(
                config["output"]["coclassify"],
                "bins_hmq/{assembler_co}.{binner_checkm}.gtdbtk.out/gtdbtk.log")
        log:
            os.path.join(config["output"]["coclassify"],
                         "logs/{assembler_co}.{binner_checkm}.gtdbtk.log")
        params:
            bin_suffix = "fa",
            out_dir = os.path.join(
                config["output"]["coclassify"],
                "bins_hmq/{assembler_co}.{binner_checkm}.gtdbtk.out")
        threads:
            config["params"]["classify"]["threads"]
        shell:
            '''
            gtdbtk classify_wf \
            --genome_dir {input}/ \
            --out_dir {params.out_dir} \
            --extension {params.bin_suffix} \
            --cpus {threads} \
            > {log}
            '''


    rule coclassify_hmq_bins_gtdbtk_all:
        input:
            expand(
                os.path.join(
                    config["output"]["coclassify"],
                    "bins_hmq/{assembler_co}.{binner_checkm}.gtdbtk.out/gtdbtk.log"),
                assembler_co=ASSEMBLERS_CO,
                binner_checkm=BINNERS_CHECKM),

            rules.cocheckm_all.input,

else:
    rule coclassify_hmq_bins_gtdbtk_all:
        input:


rule coclassify_all:
    input:
        rules.coclassify_hmq_bins_gtdbtk_all.input


rule classify_hmq_bins_gtdbtk_all:
    input:
        rules.single_classify_hmq_bins_gtdbtk_all.input,
        rules.coclassify_hmq_bins_gtdbtk_all.input


rule classify_all:
    input:
        rules.single_classify_all.input,
        rules.coclassify_all.input
