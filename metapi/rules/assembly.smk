def get_reads_for_assembly(wildcards, step, have_single=False, have_long=False):
    samples_id_list = metapi.get_samples_id_by_assembly_group(SAMPLES, wildcards.assembly_group)
    short_reads = get_short_reads_list(step, samples_id_list)

    if have_long:
        long_reads = expand(os.path.join(
            config["output"]["raw"],
            "long_reads/{sample}/{sample}.{step}{read}.fq"),
            step="raw",
            read=long_reads_suffix(),
            sample=samples_id_list)
        return [short_reads, long_reads]
    else:
        return short_reads


def assembly_input_with_short_reads(wildcards):
    if RMHOST_DO:
        return get_reads_for_assembly(wildcards, "rmhost", False, False)
    elif TRIMMING_DO:
        return get_reads_for_assembly(wildcards, "trimming", False, False)
    else:
        return get_reads_for_assembly(wildcards, "raw", False, False)


def assembly_input_with_short_and_long_reads(wildcards):
    if RMHOST_DO:
        return get_reads_for_assembly(wildcards, "rmhost", False, True)
    elif TRIMMING_DO:
        return get_reads_for_assembly(wildcards, "trimming", False, True)
    else:
        return get_reads_for_assembly(wildcards, "raw", False, True)


def get_megahit_input_str(wildcards):
    input_list = assembly_input_with_short_reads(wildcards)
    if IS_PE:
        inputstr = f'''-1 {",".join(input_list[0])} -2 {",".join(input_list[1])}'''
    else:
        inputstr = f'''-r {",".join(input_list[0])}'''
    return inputstr


if "megahit" in ASSEMBLERS:
    rule assembly_megahit:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.gfa.gz"))
        conda:
            config["envs"]["megahit"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/megahit/{assembly_group}.megahit.benchmark.txt")
        priority:
            20
        params:
            output_prefix = "{assembly_group}",
            min_contig = config["params"]["assembly"]["megahit"]["min_contig"],
            k_list = ",".join(config["params"]["assembly"]["megahit"]["k_list"]),
            presets = config["params"]["assembly"]["megahit"]["presets"],
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.megahit.out"),
            contigs = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{assembly_group}.megahit.out/{assembly_group}.contigs.fa"),
            fastg = os.path.join(config["output"]["assembly"],
                                 "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.fastg"),
            gfa = os.path.join(config["output"]["assembly"],
                               "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.gfa"),
            only_save_scaftigs = "yes" if config["params"]["assembly"]["megahit"]["only_save_scaftigs"] else "no",
            inputstr = lambda wildcards: get_megahit_input_str(wildcards)
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.megahit.log")
        shell:
            """
            set +e

            if [ -e {params.output_dir}/options.json ];
            then
                megahit --continue --out-dir {params.output_dir}
            else 
                kmeropts=""

                if [ "{params.presets}" != "" ];
                then
                    kmeropts="--presets {params.presets}"
                else
                    kmeropts="--k-list {params.k_list}"
                fi

                rm -rf {params.output_dir}

                megahit \
                {params.inputstr} \
                -t {threads} \
                $kmeropts \
                --min-contig-len {params.min_contig} \
                --out-dir {params.output_dir} \
                --out-prefix {params.output_prefix} \
                2> {log}
            fi

            knum=`grep "^>" {params.contigs} | head -1 | sed 's/>k//g' | awk -F_ '{{print $1}}'`

            megahit_toolkit contig2fastg $knum {params.contigs} > {params.fastg}

            fastg2gfa {params.fastg} > {params.gfa}
            pigz -p {threads} {params.fastg}
            pigz -p {threads} {params.gfa}

            pigz -p {threads} {params.contigs}
            mv {params.contigs}.gz {output.scaftigs}

            if [ "{params.only_save_scaftigs}" == "yes" ];
            then
                fd -t f -E "*.gz" . {params.output_dir} -x rm -rf {{}}
                rm -rf {params.output_dir}/intermediate_contigs
            fi
            """


    rule assembly_megahit_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.megahit.out/{assembly_group}.megahit.scaftigs.gfa.gz")],
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_megahit_all:
        input:


def prepare_idbaud_input(wildcards, reads):
    input_list = assembly_input_with_short_reads(wildcards)
    cmd = ""

    if IS_PE:
        if len(input_list[0]) == 1:
            cmd = f'''seqtk mergepe {input_list[0][0]} {input_list[1][0]} | seqtk seq -A - > {reads}'''
        else:
            cmd = f'''(cat {" ".join(input_list[0])} > {reads}.1.fq.gz)'''
            cmd += f''' && (cat {" ".join(input_list[1])} > {reads}.2.fq.gz)'''
            cmd += f''' && (seqtk mergepe {reads}.1.fq.gz {reads}.2.fq.gz | seqtk seq -A - > {reads})'''
            cmd += f''' && (rm -rf {reads}.1.fq.gz {reads}.2.fq.gz)'''
    else:
        if len(input_list[0]) == 1:
            cmd = f'''seqtk seq -A {input_list[0][0]} > {reads}'''
        else:
            cmd = f'''cat {" ".join(input_list[0])} > {reads}.fq.gz'''
            cmd += f''' && (seqtk seq -A {reads}.fq.gz > {reads})'''
            cmd += f''' && (rm -rf {reads}.fq.gz)'''
    return cmd


if "idba_ud" in ASSEMBLERS:
    rule assembly_idba_ud:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.idba_ud.out/{assembly_group}.idba_ud.scaftigs.fa.gz"))
        conda:
            config["envs"]["idbaud"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/idba_ud/{assembly_group}.idba_ud.benchmark.txt")
        priority:
            20
        params:
            prefix = "{assembly_group}",
            output_dir = os.path.join(config["output"]["assembly"], "scaftigs/{assembly_group}.idba_ud.out"),
            mink = config["params"]["assembly"]["idba_ud"]["mink"],
            maxk = config["params"]["assembly"]["idba_ud"]["maxk"],
            step = config["params"]["assembly"]["idba_ud"]["step"],
            min_contig = config["params"]["assembly"]["idba_ud"]["min_contig"],
            only_save_scaftigs = "yes" if config["params"]["assembly"]["idba_ud"]["only_save_scaftigs"] else "no",
            idbaud_cmd = lambda wildcards: prepare_idbaud_input(
                wildcards, 
                os.path.join(config["output"]["assembly"], f"reads_idbaud/{wildcards.assembly_group}/reads.fa")),
            reads = os.path.join(config["output"]["assembly"], "reads_idbaud/{assembly_group}/reads.fa")
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.idba_ud.log")
        shell:
            '''
            rm -rf {params.output_dir}
            mkdir -p {params.output_dir}
            mkdir -p $(dirname {params.reads})

            {params.idbaud_cmd}

            idba_ud \
            -r {params.reads} \
            --mink {params.mink} \
            --maxk {params.maxk} \
            --step {params.step} \
            --min_contig {params.min_contig} \
            -o {params.output_dir} \
            --num_threads {threads} \
            --pre_correction \
            > {log}

            rm -rf {params.reads}

            sed -i 's#^>#>{params.prefix}_#g' {params.output_dir}/scaffold.fa

            pigz -p {threads} {params.output_dir}/scaffold.fa
            mv {params.output_dir}/scaffold.fa.gz {output.scaftigs}

            if [ "{params.only_save_scaftigs}" == "yes" ];
            then
                find {params.output_dir} -type f ! -wholename "{output.scaftigs}" -delete
            fi
            '''


    rule assembly_idba_ud_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.idba_ud.out/{assembly_group}.idba_ud.scaftigs.fa.gz"),
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_idba_ud_all:
        input:


def parse_spades_params(params_file):
    import re
    import os

    if os.path.exists(params_file):
        with open(params_file, "r") as ih:
            cmd = ih.readline().strip()

            matches = re.match(r".*-k\t(.*?)\t--memory\t(\d+)\t--threads\t(\d+).*", cmd)
            if matches:
                kmers = str(matches.group(1))
                memory = str(matches.group(2))
                threads = str(matches.group(3))
                if "--only-assembler" in cmd:
                    return [kmers, memory, threads, "yes"]
                else:
                    return [kmers, memory, threads, "no"] 
            else:
                return [0, 0, 0, 0]
    else:
        return [0, 0, 0, 0]


def prepare_spades_input(wildcards, outdir):
    import os

    outdir = os.path.abspath(outdir)
    input_list = assembly_input_with_short_reads(wildcards)

    cmd_1 = ""
    cmd_2 = ""
    inputstr = ""

    if IS_PE:
        if len(input_list[0]) > 1:
            cmd_1 = f'''cat {" ".join(input_list[0])} > {outdir}/reads.1.fq.gz'''
            cmd_2 = f'''cat {" ".join(input_list[1])} > {outdir}/reads.2.fq.gz'''
        else:
            cmd_1 = f'''ln -s {os.path.abspath(input_list[0][0])} {outdir}/reads.1.fq.gz'''
            cmd_2 = f'''ln -s {os.path.abspath(input_list[1][0])} {outdir}/reads.2.fq.gz'''
        inputstr = f"-1 {outdir}/reads.1.fq.gz -2 {outdir}/reads.2.fq.gz"

    else:
        if len(input_list[0]) > 1:
            cmd_1 = f'''cat {" ".join(input_list[0])} > {outdir}/reads.fq.gz'''
            cmd_2 = "echo hello spades"
        else:
            cmd_1 = f'''ln -s {os.path.abspath(input_list[0][0])} {outdir}/reads.fq.gz'''
            cmd_2 = "echo hello spades"
        inputstr = f"-s {outdir}/reads.fq.gz"

    return [cmd_1, cmd_2, inputstr]


if "metaspades" in ASSEMBLERS:
    rule assembly_metaspades:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.gfa.gz"))
        conda:
            config["envs"]["spades"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/metaspades/{assembly_group}.metaspades.benchmark.txt")
        priority:
            20
        params:
            prefix = "{assembly_group}",
            memory = str(config["params"]["assembly"]["metaspades"]["memory"]),
            kmers = "auto" \
                if len(config["params"]["assembly"]["metaspades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["metaspades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.metaspades.out"),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["metaspades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = config["params"]["assembly"]["metaspades"]["only_save_scaftigs"],
            link_scaffolds = "yes" if config["params"]["assembly"]["metaspades"]["link_scaffolds"] else "no",
            tar_results = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.tar"),
            spades_params = lambda wildcards: parse_spades_params(
                os.path.join(os.path.join(config["output"]["assembly"],
                             "scaftigs/{wildcards.assembly_group}.metaspades.out"),
                             "params.txt")),
            spades_cmd = lambda wildcards: prepare_spades_input(
                wildcards, 
                os.path.join(config["output"]["assembly"],
                             f"reads_metaspades/{wildcards.assembly_group}")),
            reads_dir = os.path.join(config["output"]["assembly"], "reads_metaspades/{assembly_group}"),
            pe = "pe" if IS_PE else "se"
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.metaspades.log")
        shell:
            """
            if [ "{params.pe}" == "pe" ];
            then
                errorcorrect="no"
                if [ "{params.only_assembler}" != "" ];
                then
                    errorcorrect="yes"
                fi

                if [ "{params.kmers}" == "{params.spades_params[0]}" ] && [ "{params.memory}" == "{params.spades_params[1]}" ] && [ "{threads}" == "{params.spades_params[2]}" ] && [ $errorcorrect == "{params.spades_params[3]}" ];
                then
                    metaspades.py \
                    --continue \
                    -o {params.output_dir} \
                    > {log}
                else
                    rm -rf {params.output_dir}
                    mkdir -p {params.reads_dir}
                    rm -rf {params.reads_dir}/reads*

                    {params.spades_cmd[0]}
                    {params.spades_cmd[1]}

                    metaspades.py \
                    {params.spades_cmd[2]} \
                    -k {params.kmers} \
                    {params.only_assembler} \
                    --memory {params.memory} \
                    --threads {threads} \
                    --checkpoints last \
                    -o {params.output_dir} \
                    > {log}
                fi

                rm -rf {params.reads_dir}/reads*

                pigz -p {threads} {params.output_dir}/scaffolds.fasta
                mv {params.output_dir}/scaffolds.fasta.gz {params.output_dir}/{params.prefix}.metaspades.scaffolds.fa.gz

                pigz -p {threads} {params.output_dir}/contigs.fasta
                mv {params.output_dir}/contigs.fasta.gz {params.output_dir}/{params.prefix}.metaspades.contigs.fa.gz
            
                pigz -p {threads} {params.output_dir}/contigs.paths
                mv {params.output_dir}/contigs.paths.gz {params.output_dir}/{params.prefix}.metaspades.contigs.paths.gz

                pigz -p {threads} {params.output_dir}/scaffolds.paths
                mv {params.output_dir}/scaffolds.paths.gz {params.output_dir}/{params.prefix}.metaspades.scaffolds.paths.gz

                pigz -p {threads} {params.output_dir}/assembly_graph_with_scaffolds.gfa
                mv {params.output_dir}/assembly_graph_with_scaffolds.gfa.gz {output.gfa}

                if [ "{params.link_scaffolds}" == "yes" ];
                then
                    pushd {params.output_dir}
                    ln -s {params.prefix}.metaspades.scaffolds.fa.gz {params.prefix}.metaspades.scaftigs.fa.gz
                    ln -s {params.prefix}.metaspades.scaffolds.paths.gz {params.prefix}.metaspades.scaftigs.paths.gz
                    popd
                else
                    pushd {params.output_dir}
                    ln -s {params.prefix}.metaspades.contigs.fa.gz {params.prefix}.metaspades.scaftigs.fa.gz
                    ln -s {params.prefix}.metaspades.contigs.paths.gz {params.prefix}.metaspades.scaftigs.paths.gz
                    popd
                fi

                if [ "{params.only_save_scaftigs}" == "True" ];
                then
                    fd -d 1 -E "*.gz" . {params.output_dir} -x rm -rf {{}}
                else
                    rm -rf {params.output_dir}/{{corrected, misc, pipeline_state, tmp}}
                    tar -czvf {params.tar_results}.gz {params.output_dir}/K*
                fi
            else
                echo "Do not support single-end reads assembly using MetaSPAdes, you can try SPAdes or megahit, idba_ud"
                exit 1
            fi
            """


    rule assembly_metaspades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.metaspades.out/{assembly_group}.metaspades.scaftigs.gfa.gz")],
                    assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_metaspades_all:
        input:


if "spades" in ASSEMBLERS:
    rule assembly_spades:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.fa.gz")),
            gfa = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.gfa.gz"))
        conda:
            config["envs"]["spades"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/spades/{assembly_group}.spades.benchmark.txt")
        priority:
            20
        params:
            prefix = "{assembly_group}",
            memory = config["params"]["assembly"]["spades"]["memory"],
            kmers = "auto" \
                if len(config["params"]["assembly"]["spades"]["kmers"]) == 0 \
                   else ",".join(config["params"]["assembly"]["spades"]["kmers"]),
            output_dir = os.path.join(config["output"]["assembly"],
                                      "scaftigs/{assembly_group}.spades.out"),
            only_assembler = "--only-assembler" \
                if config["params"]["assembly"]["spades"]["only_assembler"] \
                   else "",
            only_save_scaftigs = config["params"]["assembly"]["spades"]["only_save_scaftigs"],
            link_scaffolds = "yes" if config["params"]["assembly"]["spades"]["link_scaffolds"] else "no",
            tar_results = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.tar"),
            spades_params = lambda wildcards: parse_spades_params(
                os.path.join(os.path.join(config["output"]["assembly"],
                             f"scaftigs/{wildcards.assembly_group}.spades.out"),
                             "params.txt")),
            spades_cmd = lambda wildcards: prepare_spades_input(
                wildcards, 
                os.path.join(config["output"]["assembly"],
                             f"reads_spades/{wildcards.assembly_group}")),
            reads_dir = os.path.join(config["output"]["assembly"], "reads_spades/{assembly_group}"),
            pe = "pe" if IS_PE else "se"
        threads:
            config["params"]["assembly"]["threads"]
        log:
            os.path.join(config["output"]["assembly"], "logs/{assembly_group}.spades.log")
        shell:
            """
            errorcorrect="no"
            if [ "{params.only_assembler}" != "" ];
            then
                errorcorrect="yes"
            fi

            if [ "{params.kmers}" == "{params.spades_params[0]}" ] && [ "{params.memory}" == "{params.spades_params[1]}" ] && [ "{threads}" == "{params.spades_params[2]}" ] && [ $errorcorrect == "{params.spades_params[3]}" ];
            then
                spades.py \
                --continue \
                -o {params.output_dir} \
                > {log}
            else
                rm -rf {params.output_dir}
                mkdir -p {params.reads_dir}
                rm -rf {params.reads_dir}/reads*

                {params.spades_cmd[0]}
                {params.spades_cmd[1]}

                spades.py \
                {params.spades_cmd[2]} \
                -k {params.kmers} \
                {params.only_assembler} \
                --memory {params.memory} \
                --threads {threads} \
                --checkpoints last \
                -o {params.output_dir} \
                > {log}
            fi

            rm -rf {params.reads_dir}/reads*

            pigz -p {threads} {params.output_dir}/scaffolds.fasta
            mv {params.output_dir}/scaffolds.fasta.gz {params.output_dir}/{params.prefix}.spades.scaffolds.fa.gz

            pigz -p {threads} {params.output_dir}/contigs.fasta
            mv {params.output_dir}/contigs.fasta.gz {params.output_dir}/{params.prefix}.spades.contigs.fa.gz
        
            pigz -p {threads} {params.output_dir}/contigs.paths
            mv {params.output_dir}/contigs.paths.gz {params.output_dir}/{params.prefix}.spades.contigs.paths.gz

            pigz -p {threads} {params.output_dir}/scaffolds.paths
            mv {params.output_dir}/scaffolds.paths.gz {params.output_dir}/{params.prefix}.spades.scaffolds.paths.gz

            pigz -p {threads} {params.output_dir}/assembly_graph_with_scaffolds.gfa
            mv {params.output_dir}/assembly_graph_with_scaffolds.gfa.gz {output.gfa}

            if [ "{params.link_scaffolds}" == "yes" ];
            then
                pushd {params.output_dir}
                ln -s {params.prefix}.spades.scaffolds.fa.gz {params.prefix}.spades.scaftigs.fa.gz
                ln -s {params.prefix}.spades.scaffolds.paths.gz {params.prefix}.spades.scaftigs.paths.gz
                popd
            else
                pushd {params.output_dir}
                ln -s {params.prefix}.spades.contigs.fa.gz {params.prefix}.spades.scaftigs.fa.gz
                ln -s {params.prefix}.spades.contigs.paths.gz {params.prefix}.spades.scaftigs.paths.gz
                popd
            fi

            if [ "{params.only_save_scaftigs}" == "True" ];
            then
                fd -d 1 -E "*.gz" . {params.output_dir} -x rm -rf {{}}
            else
                rm -rf {params.output_dir}/{{corrected, misc, pipeline_state, tmp}}
                tar -czvf {params.tar_results}.gz {params.output_dir}/K*
            fi
            """


    rule assembly_spades_all:
        input:
            expand([
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.fa.gz"),
                os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.spades.out/{assembly_group}.spades.scaftigs.gfa.gz")],
                    assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_spades_all:
        input:


if "plass" in ASSEMBLERS:
    rule assembly_plass:
        input:
            unpack(assembly_input_with_short_reads)
        output:
            proteins = os.path.join(
                config["output"]["assembly"],
                "proteins/{assembly_group}.plass.out/{assembly_group}.plass.proteins.fa.gz"),
            tmp = directory(temp(os.path.join(
                config["output"]["assembly"],
                "proteins/{assembly_group}.plass.out.tmp")))
        conda:
            config["envs"]["plass"]
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/plass/{assembly_group}.plass.benchmark.txt")
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.plass.log")
        threads:
            config["params"]["assembly"]["threads"]
        params:
            min_seq_id = config["params"]["assembly"]["plass"]["min_seq_id"],
            min_length = config["params"]["assembly"]["plass"]["min_length"],
            evalue = config["params"]["assembly"]["plass"]["evalue"],
            filter_proteins = config["params"]["assembly"]["plass"]["filter_proteins"]
        shell:
            '''
            plass assemble \
            {input} \
            --threads {threads} \
            --compressed 1 \
            --min-seq-id {params.min_seq_id} \
            --min-length {params.min_length} \
            -e {params.evalue} \
            --filter-proteins {params.filter_proteins} \
            {output.tmp}
            '''


    rule assembly_plass_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "proteins/{assembly_group}.plass.out/{assembly_group}.plass.proteins.fa.gz"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_plass_all:
        input:


def opera_ms_scaftigs_input(wildcards):
    return expand(
        os.path.join(
            config["output"]["assembly"],
            "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.megahit.scaftigs.fa.gz"),
        sample=wildcards.sample,
        assembler=config["params"]["assembly"]["opera_ms"]["short_read_assembler"])


if "opera_ms" in ASSEMBLERS:
    rule assembly_opera_ms:
        input:
            reads = assembly_input_with_short_and_long_reads,
            scaftigs = opera_ms_scaftigs_input
        output:
            scaftigs = protected(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.opera_ms.out/{assembly_group}.opera_ms.scaftigs.fa.gz"))
        benchmark:
            os.path.join(config["output"]["assembly"],
                         "benchmark/opera_ms/{assembly_group}.opera_ms.benchmark.txt")
        log:
            os.path.join(config["output"]["assembly"],
                         "logs/{assembly_group}.opera_ms.log")
        params:
            opera_ms = config["params"]["assembly"]["opera_ms"]["path"],
            prefix = "{assembly_group}",
            out_dir = os.path.join(config["output"]["assembly"],
                                   "scaftigs/{assembly_group}.opera_ms.out"),
            no_ref_clustering = "--no-ref-clustering" \
                if config["params"]["assembly"]["opera_ms"]["no_ref_clustering"] else "",
            no_strain_clustering = "--no-strain-clustering" \
                if config["params"]["assembly"]["opera_ms"]["no_strain_clustering"] else "",
            no_gap_filling = "--no-gap-filling" \
                if config["params"]["assembly"]["opera_ms"]["no_gap_filling"] else "",
            polishing = "--polishing" \
                if config["params"]["assembly"]["opera_ms"]["polishing"] else "",
            long_read_mapper = config["params"]["assembly"]["opera_ms"]["long_read_mapper"],
            short_read_assembler = config["params"]["assembly"]["opera_ms"]["short_read_assembler"],
            contig_len_threshold = config["params"]["assembly"]["opera_ms"]["contig_len_threshold"],
            contig_edge_len = config["params"]["assembly"]["opera_ms"]["contig_edge_len"],
            contig_window_len = config["params"]["assembly"]["opera_ms"]["contig_window_len"],
            genome_db = "--genome-db %s" % config["params"]["assembly"]["opera_ms"]["genome_db"] \
                if not config["params"]["assembly"]["opera_ms"]["no_ref_clustering"] \
                   else ""
        threads:
            config["params"]["assembly"]["threads"]
        run:
            shell(
                '''
                rm -rf {params.out_dir}
                pigz -p {threads} -dc {input.scaftigs} > {input.scaftigs}.fa

                perl {params.opera_ms} \
                --short-read1 {input.reads[0]} \
                --short-read2 {input.reads[1]} \
                --long-read {input.reads[2]} \
                --contig-file {input.scaftigs}.fa \
                --num-processors {threads} \
                --out-dir {params.out_dir} \
                {params.no_ref_clustering} \
                {params.no_strain_clustering} \
                {params.no_gap_filling} \
                {params.polishing} \
                {params.genome_db} \
                --long-read-mapper {params.long_read_mapper} \
                --short-read-assembler {params.short_read_assembler} \
                --contig-len-thr {params.contig_len_threshold} \
                --contig-edge-len {params.contig_edge_len} \
                --contig-window-len {params.contig_window_len} \
                >{log} 2>&1
                ''')

            shell('''pigz -p {threads} {params.out_dir}/contigs.fasta''')

            if params.polishing != "":
                shell(
                    '''
                    pigz -p {threads} {params.out_dir}/contigs.polished.fasta

                    pushd {params.out_dir} && \
                    ln -s contigs.polished.fasta.gz {params.prefix}.opera_ms.scaftigs.fa.gz && \
                    popd
                    ''')
            else:
                shell(
                    '''
                    pushd {params.out_dir} && \
                    ln -s contigs.fasta.gz {params.prefix}.opera_ms.scaftigs.fa.gz && \
                    popd
                    ''')
            shell('''rm -rf {input.scaftigs}.fa''')


    rule assembly_opera_ms_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.opera_ms.out/{assembly_group}.opera_ms.scaftigs.fa.gz"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

else:
    rule assembly_opera_ms_all:
        input:


if len(ASSEMBLERS) != 0:
    if config["params"]["assembly"]["metaquast"]["do"] and IS_PE:
        rule assembly_metaquast:
            input:
                reads = assembly_input_with_short_reads,
                scaftigs = os.path.join(
                    config["output"]["assembly"],
                    "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
            output:
                protected(os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{assembly_group}.{assembler}.metaquast.out/combined_reference/report.tsv"))
            conda:
                config["envs"]["quast"]
            benchmark:
                os.path.join(config["output"]["assembly"],
                             "benchmark/metaquast_{assembler}/{assembly_group}.{assembler}.metaquast.benchmark.txt")

            log:
                os.path.join(config["output"]["assembly"],
                             "logs/{assembly_group}.{assembler}.metaquast.log")
            params:
                labels = "{assembly_group}.{assembler}",
                output_dir = os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{assembly_group}.{assembler}.metaquast.out")
            threads:
                config["params"]["assembly"]["metaquast"]["threads"]
            shell:
                '''
                metaquast.py \
                {input.scaftigs} \
                --pe1 {input.reads[0]} \
                --pe2 {input.reads[1]} \
                --output-dir {params.output_dir} \
                --labels {params.labels} \
                --circos \
                --ran-finding \
                --conserved-genes-finding \
                --threads {threads} \
                2> {log}
                '''


        rule assembly_metaquast_multiqc:
            input:
                expand(os.path.join(
                    config["output"]["assembly"],
                    "metaquast/{assembly_group}.{{assembler}}.metaquast.out/combined_reference/report.tsv"),
                    assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
            output:
                html = os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                data_dir = directory(
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report_data"))
            conda:
                config["envs"]["multiqc"]
            log:
                os.path.join(config["output"]["assembly"], "logs/multiqc_{assembler}_metaquast.log")
            params:
                output_dir = os.path.join(
                    config["output"]["assembly"],
                    "report/{assembler}_metaquast")
            shell:
                '''
                multiqc \
                --outdir {params.output_dir} \
                --title metaquast \
                --module quast \
                {input} \
                2> {log}
                '''


        rule assembly_metaquast_all:
            input:
                expand([
                    os.path.join(
                        config["output"]["assembly"],
                        "metaquast/{assembly_group}.{assembler}.metaquast.out/combined_reference/report.tsv"),
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report.html"),
                    os.path.join(
                        config["output"]["assembly"],
                        "report/{assembler}_metaquast/metaquast_multiqc_report_data")],
                       assembler=ASSEMBLERS,
                       assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

    else:
        rule assembly_metaquast_all:
            input:

           
    rule assembly_report:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz")
        output:
            report = os.path.join(
                config["output"]["assembly"],
                "report/{assembler}_stats/{assembly_group}.{assembler}.scaftigs.seqtk.comp.tsv.gz")
        conda:
            config["envs"]["report"]
        priority:
            25
        params:
            assembly_group = "{assembly_group}",
            assembler = "{assembler}"
        shell:
            '''
            seqtk comp {input.scaftigs} | \
            awk \
            'BEGIN \
            {{print "assembly_group\tassembler\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts"}}; \
            {{print "{params.assembly_group}" "\t" "{params.assembler}" "\t" $0}}' | \
            gzip -c > {output.report}
            '''


    rule assembly_report_merge:
        input:
            comp_list = expand(
                os.path.join(
                    config["output"]["assembly"],
                    "report/{{assembler}}_stats/{assembly_group}.{{assembler}}.scaftigs.seqtk.comp.tsv.gz"),
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)
        output:
            summary = os.path.join(
                config["output"]["assembly"],
                "report/assembly_stats_{assembler}.tsv")
        params:
            min_length = config["params"]["assembly"]["report"]["min_length"],
            len_ranges = config["params"]["assembly"]["report"]["len_ranges"]
        threads:
            config["params"]["assembly"]["threads"]
        run:
            comp_list = [(i, params.min_length) for i in input.comp_list]
            metapi.assembler_init(params.len_ranges, ["assembly_group", "assembler"])
            metapi.merge(comp_list, metapi.parse_assembly,
                         threads, output=output.summary)


    rule assembly_report_all:
        input:
            expand(os.path.join(
                config["output"]["assembly"],
                "report/assembly_stats_{assembler}.tsv"),
                   assembler=ASSEMBLERS)

else:
    rule assembly_report_all:
        input:


    rule assembly_metaquast_all:
        input:


rule assembly_all:
    input:
        rules.assembly_megahit_all.input,
        rules.assembly_idba_ud_all.input,
        rules.assembly_metaspades_all.input,
        rules.assembly_spades_all.input,
        rules.assembly_plass_all.input,
        rules.assembly_opera_ms_all.input,

        rules.assembly_report_all.input,
        rules.assembly_metaquast_all.input#,

        #rules.rmhost_all.input,
        #rules.qcreport_all.input


localrules:
    assembly_spades_all,
    assembly_idba_ud_all,
    assembly_megahit_all,
    assembly_metaspades_all,
    assembly_opera_ms_all,
    assembly_plass_all,
    assembly_metaquast_all,
    assembly_report_all,
    assembly_all