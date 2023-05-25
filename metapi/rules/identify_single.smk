rule identify_virsorter2_setup_db:
    output:
        directory(os.path.join(config["params"]["identify"]["virsorter2"]["db"], "hmm")),
        directory(os.path.join(config["params"]["identify"]["virsorter2"]["db"], "rbs")),
        directory(os.path.join(config["params"]["identify"]["virsorter2"]["db"], "group"))
    log:
        os.path.join(config["output"]["identify"], "logs/identify_virsorter2_setup_db/identify_virsorter2_setup_db.log")
    benchmark:
        os.path.join(config["output"]["identify"], "benchmark/identify_virsorter2_setup_db/identify_virsorter2_setup_db.txt")
    params:
        db_dir = config["params"]["identify"]["virsorter2"]["db"]
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["virsorter2"]
    shell:
        '''
        mkdir -p {params.db_dir}

        virsorter setup --db-dir {params.db_dir} --jobs {threads} >{log} 2>&1
        '''


# https://github.com/EddyRivasLab/hmmer/issues/161
# hmmsearch threads: 2 (recommand)
rule identify_virsorter2_config:
    input:
        expand(os.path.join(
            config["params"]["identify"]["virsorter2"]["db"], "{dbdir}"),
            dbdir=["hmm", "rbs", "group"])
    output:
        os.path.join(config["output"]["identify"], "config/virsorter2-template-config.yaml")
    log:
        os.path.join(config["output"]["identify"], "logs/identify_virsorter2_config/identify_virsorter2_config.log")
    benchmark:
        os.path.join(config["output"]["identify"], "benchmark/identify_virsorter2_config/identify_virsorter2_config.txt")
    params:
        db_dir = config["params"]["identify"]["virsorter2"]["db"]
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["virsorter2"]
    shell:
        '''
        configfile=`python -c 'import os,virsorter;print(os.path.join(virsorter.__path__[0], "template-config.yaml"))'`

        if [ -f $configfile ];
        then
            cp $configfile {output}
        else
            virsorter config --init-source --db-dir={params.db_dir} >{log} 2>&1

            virsorter config --set GENERAL_THREADS={threads} >>{log} 2>&1
            virsorter config --set HMMSEARCH_THREADS=2 >>{log} 2>&1
            virsorter config --set CLASSIFY_THREADS={threads} >>{log} 2>&1

            cp $configfile {output}
        fi
        '''


checkpoint identify_virsorter2_prepare:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        scaftigs_dir = directory(os.path.join(
            config["output"]["assembly"],
            "scaftigs_splited/{binning_group}.{assembly_group}.{assembler}"))
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_virsorter2_prepare/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_virsorter2_prepare/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        split_contigs_num = config["params"]["identify"]["virsorter2"]["split_contigs_num"]
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["report"]
    shell:
        '''
        rm -rf {output.scaftigs_dir}

        seqkit split {input.scaftigs} \
        --by-size {params.split_contigs_num} \
        --line-width 0 \
        --out-dir {output.scaftigs_dir} \
        --threads {threads} \
        >{log} 2>&1
        '''


# Avoid creating the conda_envs directory at the same time when run virsorter2 using split scaftigs
rule identify_virsorter2_init_run:
    input:
        config_file = os.path.join(config["output"]["identify"], "config/virsorter2-template-config.yaml"),
        lambda_virus = os.path.join(DATA_DIR, "lambda_virus.fa")
    output:
        init_success = os.path.join(config["output"]["identify"], "config/virsorter2-init-run-success")
    log:
        os.path.join(config["output"]["identify"], "logs/identify_virsorter2_init_run/identify_virsorter2_init_run.log")
    benchmark:
        os.path.join(config["output"]["identify"], "benchmark/identify_virsorter2_init_run/identify_virsorter2_init_run.txt")
    params:
        label = "lambda_virus",
        working_dir = os.path.join(config["output"]["identify"], "config/lambda_virus.vs2.out"),
        include_groups = ",".join(config["params"]["identify"]["virsorter2"]["include_groups"]),
        min_length = config["params"]["identify"]["virsorter2"]["min_length"],
        min_score = config["params"]["identify"]["virsorter2"]["min_score"],
        provirus_off = "--provirus-off" if config["params"]["identify"]["virsorter2"]["provirus_off"] else "",
        prep_for_dramv = "--prep-for-dramv" if config["params"]["identify"]["virsorter2"]["prep_for_dramv"] else "",
        rm_tmpdir = "--rm-tmpdir" if config["params"]["identify"]["virsorter2"]["rm_tmpdir"] else "",
        keep_original_seq = "--keep-original-seq" if config["params"]["identify"]["virsorter2"]["keep_original_seq"] else ""
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["virsorter2"]
    shell:
        '''
        rm -rf {params.working_dir}
        mkdir -p {params.working_dir}

        virsorter run \
        {params.prep_for_dramv} \
        {params.provirus_off} \
        {params.rm_tmpdir} \
        {params.keep_original_seq} \
        --working-dir {params.working_dir} \
        --seqfile {input.lambda_virus} \
        --label {params.label} \
        --include-groups {params.include_groups} \
        --min-length {params.min_length} \
        --min-score {params.min_score} \
        --jobs {threads} all \
        >{log} 2>&1

        touch {output}
        '''


localrules:
    identify_virsorter2_setup_db,
    identify_virsorter2_config,
    identify_virsorter2_init_run


rule identify_virsorter2:
    input:
        init_success = os.path.join(config["output"]["identify"], "config/virsorter2-init-run-success"),
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs_splited/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.part_{split_num}.fa.gz")
    output:
        os.path.join(
            config["output"]["identify"],
            "vmags/{binning_group}.{assembly_group}.{assembler}/virsorter2/virsorter2_{split_num}/virsorter2_done")
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_virsorter2/{binning_group}.{assembly_group}.{assembler}.{split_num}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_virsorter2/{binning_group}.{assembly_group}.{assembler}.{split_num}.txt")
    params:
        label = "{binning_group}.{assembly_group}.{assembler}",
        include_groups = ",".join(config["params"]["identify"]["virsorter2"]["include_groups"]),
        working_dir = os.path.join(config["output"]["identify"], "vmags/{binning_group}.{assembly_group}.{assembler}/virsorter2/virsorter2_{split_num}"),
        min_length = config["params"]["identify"]["virsorter2"]["min_length"],
        min_score = config["params"]["identify"]["virsorter2"]["min_score"],
        provirus_off = "--provirus-off" if config["params"]["identify"]["virsorter2"]["provirus_off"] else "",
        prep_for_dramv = "--prep-for-dramv" if config["params"]["identify"]["virsorter2"]["prep_for_dramv"] else "",
        rm_tmpdir = "--rm-tmpdir" if config["params"]["identify"]["virsorter2"]["rm_tmpdir"] else "",
        keep_original_seq = "--keep-original-seq" if config["params"]["identify"]["virsorter2"]["keep_original_seq"] else ""
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["virsorter2"]
    shell:
        '''
        rm -rf {params.working_dir}
        mkdir -p {params.working_dir}

        set +e

        virsorter run \
        {params.prep_for_dramv} \
        {params.provirus_off} \
        {params.rm_tmpdir} \
        {params.keep_original_seq} \
        --working-dir {params.working_dir} \
        --seqfile {input.scaftigs} \
        --label {params.label} \
        --include-groups {params.include_groups} \
        --min-length {params.min_length} \
        --min-score {params.min_score} \
        --jobs {threads} all \
        >{log} 2>&1

        exitcode=$?
        echo "Exit code is: $exitcode" >> {log}

        if [ $exitcode -eq 1 ];
        then
            grep -oEi "No genes from the contigs are left in iter-0/all.pdg.faa after preprocess" {log} 
            grepcode=$?
            if [ $grepcode -eq 0 ];
            then
                echo "Touch {output}" >> {log}
                touch {output} >> {log} 2>&1
                exit 0
            else
                grep -oEi "Error in rule circular_linear_split" {log}
                grepcode=$?
                if [ $grepcode -eq 0 ];
                then
                    echo "Touch {output}" >> {log}
                    touch {output} >> {log} 2>&1
                    exit 0
                else
                    echo "Runing failed, check Virsorter log please." >> {log} 2>&1
                    exit $exitcode
                fi
            fi
        else
            echo "Touch {output}" >> {log}
            touch {output} >> {log} 2>&1
        fi
        '''


def aggregate_identify_virsorter2_output(wildcards):
    checkpoint_output = checkpoints.identify_virsorter2_prepare.get(**wildcards).output.scaftigs_dir
    label = wildcards.binning_group + "." + wildcards.assembly_group + "." + wildcards.assembler

    return expand(os.path.join(
        config["output"]["identify"],
        "vmags/{binning_group}.{assembly_group}.{assembler}/virsorter2/virsorter2_{split_num}/virsorter2_done"),
        binning_group=wildcards.binning_group,
        assembly_group=wildcards.assembly_group,
        assembler=wildcards.assembler,
        split_num=list(set([i.split("/")[0] \
            for i in glob_wildcards(os.path.join(
                checkpoint_output,
                f"{label}.scaftigs.part_{{split_num}}.fa.gz")).split_num])))


rule identify_virsorter2_merge:
    input:
        virsorter2_done = aggregate_identify_virsorter2_output
    output:
        expand(os.path.join(
            config["output"]["identify"],
            "vmags/{{binning_group}}.{{assembly_group}}.{{assembler}}/virsorter2/{{binning_group}}.{{assembly_group}}.{{assembler}}.virsorter2.{suffix}.gz"),
            suffix=["combined.fa", "score.tsv", "boundary.tsv"])
    params:
        label = "{binning_group}.{assembly_group}.{assembler}"
    threads:
        config["params"]["identify"]["threads"]
    run:
        import os
        import subprocess

        combined_fa_list = []
        score_tsv_list = []
        boundary_tsv_list = []

        for i in input.virsorter2_done:
            vs2_dir = os.path.dirname(i)
            combined_fa = os.path.join(vs2_dir, params.label + "-final-viral-combined.fa")
            score_tsv = os.path.join(vs2_dir, params.label + "-final-viral-score.tsv")
            boundary_tsv = os.path.join(vs2_dir, params.label + "-final-viral-boundary.tsv")

            if os.path.exists(combined_fa) and (os.path.getsize(combined_fa) > 0):
                combined_fa_list.append(combined_fa)
                if os.path.exists(score_tsv):
                    score_tsv_list.append(score_tsv)
                if os.path.exists(boundary_tsv):
                    boundary_tsv_list.append(boundary_tsv)

        if len(combined_fa_list) > 0:
            fa_str = " ".join(combined_fa_list)
            subprocess.run(f"cat {fa_str} | pigz -cf > {output[0]}", shell=True)
        else:
            subprocess.run(f"touch {output[0]}", shell=True)

        if len(score_tsv_list) > 0:
            metapi.merge(score_tsv_list, metapi.parse, threads, output=output[1])
        else:
            subprocess.run(f"touch {output[1]}", shell=True)

        if len(boundary_tsv_list) > 0:
            metapi.merge(boundary_tsv_list, metapi.parse, threads, output=output[2])
        else:
            subprocess.run(f"touch {output[2]}", shell=True)


if config["params"]["identify"]["virsorter2"]["do"]:
    rule identify_virsorter2_all:
        input:
            expand(expand(
                os.path.join(config["output"]["identify"],
                "vmags/{binning_group}.{assembly_group}.{assembler}/virsorter2/{binning_group}.{assembly_group}.{assembler}.virsorter2.{{suffix}}.gz"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
                suffix=["combined.fa", "score.tsv", "boundary.tsv"])

else:
    rule identify_virsorter2_all:
        input:


rule identify_deepvirfinder:
    input:
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz")
    output:
        expand(os.path.join(
            config["output"]["identify"],
            "vmags/{{binning_group}}.{{assembly_group}}.{{assembler}}/deepvirfinder/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz_gt{min_length}bp_dvfpred.txt.gz"),
            min_length=config["params"]["identify"]["deepvirfinder"]["min_length"])
    log:
        os.path.join(
            config["output"]["identify"],
            "logs/identify_deepvirfinder/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["identify"],
            "benchmark/identify_deepvirfinder/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        deepvirfinder = config["params"]["identify"]["deepvirfinder"]["script"],
        out_dir = os.path.join(config["output"]["identify"], "vmags/{binning_group}.{assembly_group}.{assembler}/deepvirfinder"),
        min_length=config["params"]["identify"]["deepvirfinder"]["min_length"]
    threads:
        config["params"]["identify"]["threads"]
    conda:
        config["envs"]["deepvirfinder"]
    shell:
        '''
        rm -rf {params.out_dir}

        set +e

        python {params.deepvirfinder} \
        --in {input} \
        --out {params.out_dir} \
        --len {params.min_length} \
        --core {threads} \
        >{log} 2>&1

        exitcode=$?
        echo "Exit code is: $exitcode" >> {log} 2>&1

        DVFGZ={output}
        DVF=${{DVFGZ%.gz}}
        if [ -f $DVF ];
        then
            pigz -f $DVF
        fi

        if [ $exitcode -eq 1 ];
        then
            grep -oEi "ValueError: not enough values to unpack" {log}
            grepcode=$?
            if [ $grepcode -eq 0 ];
            then
                touch {output} >> {log} 2>&1
                echo "Touch empty file: {output}" >> {log} 2>&1
                exit 0
            else
                echo "Runing failed, check DeepVirFinder log please." >> {log} 2>&1
                exit $exitcode
            fi
        fi
        '''


rule identify_phamb_deepvirfinder_extract_contigs:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.merged.scaftigs.fa.gz"),
        metadata = os.path.join(
            config["output"]["assembly"],
            "scaftigs_merged/{binning_group}.{assembler}/{binning_group}.{assembler}.metadata.tsv.gz"),
        dvf_anno = expand(os.path.join(
            config["output"]["identify"],
            "vmags/{{binning_group}}.{{assembly_group}}.{{assembler}}/deepvirfinder/{{binning_group}}.{{assembly_group}}.{{assembler}}.scaftigs.fa.gz_gt{min_length}bp_dvfpred.txt.gz"),
            min_length=config["params"]["identify"]["deepvirfinder"]["min_length"])
    output:
        fna = os.path.join(
            config["output"]["identify"],
            "vmags/{binning_group}.{assembly_group}.{assembler}/deepvirfinder/{binning_group}.{assembly_group}.{assembler}.deepvirfinder.combined.fa.gz")
    params:
        assembly_group = "{assembly_group}"
    run:
        from Bio import SeqIO
        import gzip

        ### record assembly_group : alias ###
        tab = pd.read_table(input.metadata)
        # vamb_id_assembly = {vamb_id : binning_assembly.split(".")[-1] for vamb_id, binning_assembly in zip(tab.iloc[:,1], tab.iloc[:, 0])}
        assembly_vamb_id = {binning_assembly.split(".")[-1] : vamb_id for vamb_id, binning_assembly in zip(tab.iloc[:,1], tab.iloc[:, 0])}
        # seems like metapi.get_samples_id_by_binning_group(SAMPLES, wildcards.binning_group) will do the trick

        # shell("zcat {input.dvf_anno} | perl -lane 's/\s.*?\t/\t/; print' > {params.temp_metadata}")
        dvf_dict = dict() #defaultdict(int)
        vamb_id = assembly_vamb_id[params.assembly_group]

        with gzip.open(input.dvf_anno[0], "rt") as records:
            records.readline()# title line
            for record in records:
                name, length, score, p = record.split("\t")

                if float(score) > 0.5 and float(p) < 0.05:
                    name = name.split(" ")[0]
                    name = vamb_id + 'C' + name
                    dvf_dict[name] = 1
                # print(record, end="")

        dvf_viral_list = list(dvf_dict.keys())
        with gzip.open(output.fna, "wt") as f:
            for record in SeqIO.parse(gzip.open(input.scaftigs, 'rt'), 'fasta'):
                if record.description in dvf_viral_list:
                    f.write(record.format("fasta"))


rule identify_phamb_deepvirfinder_merge:
    input:
        lambda wildcards: expand(os.path.join(
            config["output"]["identify"],
            "vmags/{{binning_group}}.{assembly_group}.{{assembler}}/deepvirfinder/{{binning_group}}.{assembly_group}.{{assembler}}.scaftigs.fa.gz_gt{min_length}bp_dvfpred.txt.gz"),
            min_length=config["params"]["identify"]["deepvirfinder"]["min_length"],
            assembly_group=sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, wildcards.binning_group)))
    output:
        dvf = os.path.join(config["output"]["identify"], "annotations/{binning_group}.{assembler}/all.DVF.predictions.txt.gz")
    params:
        binning_group = "{binning_group}"
    run:
        import gzip

        assembly_groups = sorted(metapi.get_assembly_group_by_binning_group(SAMPLES, params.binning_group))

        with gzip.open(output.dvf, "wt") as oh:
            assembly_group = os.path.basename(input[0]).split(".")[1]
            assembly_index = int(assembly_groups.index(assembly_group)) + 1
            assembly_group = f'''S{assembly_index}'''
            with gzip.open(input[0], "rt") as ih:
                oh.write(ih.readline())
                for line in ih:
                    oh.write(f'''{assembly_group}C{line}''')
            for dvfpred in input[1:]:
                assembly_group = os.path.basename(dvfpred).split(".")[1]
                assembly_index = int(assembly_groups.index(assembly_group)) + 1
                assembly_group = f'''S{assembly_index}'''
                with gzip.open(dvfpred, "rt") as ih:
                    ih.readline()
                    for line in ih:
                        oh.write(f'''{assembly_group}C{line}''')


if config["params"]["identify"]["deepvirfinder"]["do"]:
    rule identify_deepvirfinder_all:
        input:
            expand(os.path.join(
                config["output"]["identify"],
                "vmags/{binning_group}.{assembly_group}.{assembler}/deepvirfinder/{binning_group}.{assembly_group}.{assembler}.deepvirfinder.combined.fa.gz"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"]),
            expand(os.path.join(
                config["output"]["identify"],
                "annotations/{binning_group}.{assembler}/all.DVF.predictions.txt.gz"),
                binning_group=SAMPLES_BINNING_GROUP_LIST,
                assembler=ASSEMBLERS)

else:
    rule identify_deepvirfinder_all:
        input:


rule identify_single_all:
    input:
        rules.identify_virsorter2_all.input,
        rules.identify_deepvirfinder_all.input


localrules:
    identify_virsorter2_all,
    identify_deepvirfinder_all,
    identify_single_all
