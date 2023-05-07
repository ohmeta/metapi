rule binning_metabat2_coverage:
    input:
        bam = lambda wildcards: metapi.get_samples_bai(wildcards, SAMPLES, config["output"]["alignment"], "bam"),
        bai = lambda wildcards: metapi.get_samples_bai(wildcards, SAMPLES, config["output"]["alignment"], "bai")
    output:
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.coverage.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_metabat2_coverage/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_metabat2_coverage/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        percent_identity = config["params"]["binning"]["metabat2"]["percent_identity"],
        min_map_qual = config["params"]["binning"]["metabat2"]["min_map_qual"],
        output_paired_contigs = "--pairedContigs %s" % \
            os.path.join(
                config["output"]["binning"],
                "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.paired_contigs") \
                if config["params"]["binning"]["metabat2"]["output_paired_contigs"] \
                    else "",
        output_gc = "--outputGC %s" % \
            os.path.join(
                config["output"]["binning"],
                "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.gc") \
                if config["params"]["binning"]["metabat2"]["output_gc"] \
                    else "",
        output_gc_window = "--gcWindow %s" % \
            os.path.join(
                config["output"]["binning"],
                "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.gc_window") \
                if config["params"]["binning"]["metabat2"]["output_gc_window"] \
                    else ""
    priority:
        28
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["metabat2"]
    shell:
        '''
        COVERAGE={output.coverage}

        jgi_summarize_bam_contig_depths \
        --outputDepth ${{COVERAGE%.gz}} \
        --percentIdentity {params.percent_identity} \
        --minMapQual {params.min_map_qual} \
        {params.output_paired_contigs} \
        {params.output_gc} \
        {params.output_gc_window} \
        {input.bam} \
        >{log} 2>&1

        pigz -f ${{COVERAGE%.gz}}
        '''


rule binning_metabat2_coverage_all:
    input:
        expand(os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.coverage.gz"),
            zip,
            binning_group=ASSEMBLY_GROUPS["binning_group"],
            assembly_group=ASSEMBLY_GROUPS["assembly_group"],
            assembler=ASSEMBLY_GROUPS["assembler"])


localrules:
    binning_metabat2_coverage_all


rule binning_metabat2:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.coverage.gz")
    output:
        os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/metabat2/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_metabat2/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_metabat2/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        mags_dir = os.path.join(config["output"]["binning"], "mags/{binning_group}.{assembly_group}.{assembler}/metabat2"),
        bin_prefix = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/metabat2/{binning_group}.{assembly_group}.{assembler}.metabat2.bin"),
        min_contig = config["params"]["binning"]["metabat2"]["min_contig"],
        max_p = config["params"]["binning"]["metabat2"]["maxP"],
        min_s = config["params"]["binning"]["metabat2"]["minS"],
        max_edges = config["params"]["binning"]["metabat2"]["maxEdges"],
        p_tnf = config["params"]["binning"]["metabat2"]["pTNF"],
        no_add = "--noAdd" if config["params"]["binning"]["metabat2"]["noAdd"] else "",
        min_cv = config["params"]["binning"]["metabat2"]["minCV"],
        min_cv_sum = config["params"]["binning"]["metabat2"]["minCVSum"],
        min_cls_size = config["params"]["binning"]["metabat2"]["minClsSize"],
        save_cls = "--saveCls" \
            if config["params"]["binning"]["metabat2"]["saveCls"] else "",
        seed = config["params"]["binning"]["metabat2"]["seed"]
    priority:
        30
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["metabat2"]
    shell:
        '''
        rm -rf {params.mags_dir}

        COVERAGE={input.coverage}
        pigz -p {threads} -dkf $COVERAGE

        metabat2 \
        --inFile {input.scaftigs} \
        --abdFile ${{COVERAGE%.gz}} \
        --outFile {params.bin_prefix} \
        --minContig {params.min_contig} \
        --maxP {params.max_p} \
        --minS {params.min_s} \
        --maxEdges {params.max_edges} \
        --pTNF {params.p_tnf} \
        {params.no_add} \
        --minCV {params.min_cv} \
        --minCVSum {params.min_cv_sum} \
        {params.save_cls} \
        --seed {params.seed} \
        --numThreads {threads} \
        --verbose \
        >{log} 2>&1

        rm -rf ${{COVERAGE%.gz}}

        for FILESTR in `ls {params.bin_prefix}*`
        do
            if [ -f $FILESTR ] && [ "${{FILESTR##*.}}" != "gz" ];
            then
                pigz -f -p {threads} $FILESTR
            fi
        done

        touch {output}
        '''


if config["params"]["binning"]["metabat2"]["do"]:
    rule binning_metabat2_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/metabat2/binning_done"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"])

else:
    rule binning_metabat2_all:
        input:


rule binning_maxbin2_coverage:
    input:
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.metabat2.coverage.gz")
    output:
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.maxbin2.coverage.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_maxbin2_coverage/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_maxbin2_coverage/{binning_group}.{assembly_group}.{assembler}.txt")
    priority:
        30
    threads:
        1
    shell:
        '''
        zcat {input.coverage} | \
        cut -f1,3 | tail -n +2 | \
        pigz -cf \
        > {output.coverage} \
        2> {log}
        '''


rule binning_maxbin2:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.maxbin2.coverage.gz")
    output:
        os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/maxbin2/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_maxbin2/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_maxbin2/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        wrapper_dir = WRAPPER_DIR,
        mags_dir = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/maxbin2"),
        bin_prefix = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/maxbin2/{binning_group}.{assembly_group}.{assembler}.maxbin2.bin"),
        min_contig = config["params"]["binning"]["maxbin2"]["min_contig"],
        max_iteration = config["params"]["binning"]["maxbin2"]["max_iteration"],
        prob_threshold = config["params"]["binning"]["maxbin2"]["prob_threshold"],
        plotmarker = "-plotmarker" if config["params"]["binning"]["maxbin2"]["plotmarker"] \
            else "",
        markerset = config["params"]["binning"]["maxbin2"]["markerset"]
    priority:
        30
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["maxbin2"]
    shell:
        '''
        rm -rf {params.mags_dir}
        mkdir -p {params.mags_dir}

        COVERAGE={input.coverage}
        pigz -dkf -p {threads} $COVERAGE

        set +e
        run_MaxBin.pl \
        -thread {threads} \
        -contig {input.scaftigs} \
        -abund ${{COVERAGE%.gz}} \
        -min_contig_length {params.min_contig} \
        -max_iteration {params.max_iteration} \
        -prob_threshold {params.prob_threshold} \
        {params.plotmarker} \
        -markerset {params.markerset} \
        -out {params.bin_prefix} \
        >{log} 2>&1

        rm -rf ${{COVERAGE%.gz}}

        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            grep -oEi 'Program stop' {log}
            grepcode=$?
            if [ $grepcode -eq 0 ]
            then
                exit 0
            else
                exit $exitcode
            fi
        fi

        python {params.wrapper_dir}/maxbin2_postprocess.py \
        {params.mags_dir}

        for FILESTR in `ls {params.bin_prefix}*`
        do
            if [ -f $FILESTR ] && [ "${{FILESTR##*.}}" != "gz" ];
            then
                pigz -f -p {threads} $FILESTR
            fi
        done

        touch {output}
        '''


if config["params"]["binning"]["maxbin2"]["do"]:
    rule binning_maxbin2_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/maxbin2/binning_done"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"])

else:
    rule binning_maxbin2_all:
        input:


'''
if config["params"]["binning"]["canopy"]["do"]:
    rule binning_canopy_coverage:
        input:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}/{sample}.{assembler}.metabat2.coverage")
        output:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}/{sample}.{assembler}.canopy.coverage")
        priority:
            30
        run:
            import pandas as pd

            df = pd.read_csv(input.coverage, sep='\t')
            df.iloc[:, [0, 3]].to_csv(output.coverage, header=None, sep='\t', index=False)
'''


rule binning_concoct_cut_bed:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
    output:
        scaftigs_cut = os.path.join(
            config["output"]["assembly"],
            "scaftigs_cut/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.cut.fa.gz"),
        scaftigs_bed = os.path.join(
            config["output"]["assembly"],
            "scaftigs_cut/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.cut.bed.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_concoct_cut_bed/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_concoct_cut_bed/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        chunk_size = config["params"]["binning"]["concoct"]["chunk_size"],
        overlap_size = config["params"]["binning"]["concoct"]["overlap_size"]
    threads:
        1
    conda:
        config["envs"]["concoct"]
    shell:
        '''
        rm -rf {output.scaftigs_cut}
        rm -rf {output.scaftigs_bed}

        SCAFTIGS={input.scaftigs}
        BED={output.scaftigs_bed}

        pigz -dkf $SCAFTIGS

        cut_up_fasta.py \
        ${{SCAFTIGS%.gz}} \
        --chunk_size {params.chunk_size} \
        --overlap_size {params.overlap_size} \
        --merge_last \
        --bedfile ${{BED%.gz}} | \
        pigz -cf > {output.scaftigs_cut} \
        2>{log}

        rm -rf ${{SCAFTIGS%.gz}}
        pigz -f ${{BED%.gz}}
        '''


rule binning_concoct_coverage:
    input:
        bam = lambda wildcards: metapi.get_samples_bai(wildcards, SAMPLES, config["output"]["alignment"], "bam"),
        bai = lambda wildcards: metapi.get_samples_bai(wildcards, SAMPLES, config["output"]["alignment"], "bai"),
        scaftigs_bed = os.path.join(
            config["output"]["assembly"],
            "scaftigs_cut/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.cut.bed.gz"),
    output:
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.concoct.coverage.gz")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_concoct_coverage/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/binning_concoct_coverage/{binning_group}.{assembly_group}.{assembler}.txt")
    priority:
        30
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["concoct"]
    shell:
        '''
        BED={input.scaftigs_bed}
        pigz -dkf $BED

        concoct_coverage_table.py \
        ${{BED%.gz}} \
        {input.bam} | \
        pigz -cf > {output.coverage} \
        2> {log}

        rm -rf ${{BED%.gz}}
        '''


rule binning_concoct:
    input:
        scaftigs = os.path.join(
            config["output"]["assembly"],
            "scaftigs/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.fa.gz"),
        scaftigs_cut = os.path.join(
            config["output"]["assembly"],
            "scaftigs_cut/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.cut.fa.gz"),
        scaftigs_bed = os.path.join(
            config["output"]["assembly"],
            "scaftigs_cut/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.scaftigs.cut.bed.gz"),
        coverage = os.path.join(
            config["output"]["binning"],
            "coverage/{binning_group}.{assembly_group}.{assembler}/{binning_group}.{assembly_group}.{assembler}.concoct.coverage.gz")
    output:
        os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/concoct/binning_done")
    log:
        os.path.join(
            config["output"]["binning"],
            "logs/binning_concoct/{binning_group}.{assembly_group}.{assembler}.log")
    benchmark:
        os.path.join(
            config["output"]["binning"],
            "benchmark/concoct/{binning_group}.{assembly_group}.{assembler}.txt")
    params:
        wrapper_dir = WRAPPER_DIR,
        mags_dir = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/concoct"),
        clusters = config["params"]["binning"]["concoct"]["clusters"],
        kmer_length = config["params"]["binning"]["concoct"]["kmer_length"],
        length_threshold = config["params"]["binning"]["concoct"]["length_threshold"],
        read_length = config["params"]["binning"]["concoct"]["read_length"],
        total_percentage_pca = config["params"]["binning"]["concoct"]["total_percentage_pca"],
        iterations = config["params"]["binning"]["concoct"]["iterations"],
        seed = config["params"]["binning"]["concoct"]["seed"],
        no_cov_normalization = "--no_cov_normalization" \
            if config["params"]["binning"]["concoct"]["no_cov_normalization"] \
                else "",
        no_total_coverage = "--no_total_coverage" \
            if config["params"]["binning"]["concoct"]["no_total_coverage"] \
                else "",
        no_original_data = "--no_original_data" \
            if config["params"]["binning"]["concoct"]["no_original_data"] \
                else "",
        coverage_out = "--coverage_out" \
            if config["params"]["binning"]["concoct"]["coverage_out"] \
                else "",
        bin_prefix = os.path.join(
            config["output"]["binning"],
            "mags/{binning_group}.{assembly_group}.{assembler}/concoct/{binning_group}.{assembly_group}.{assembler}.concoct.bin"),
    threads:
        config["params"]["binning"]["threads"]
    conda:
        config["envs"]["concoct"]
    shell:
        '''
        rm -rf {params.mags_dir}
        mkdir -p {params.mags_dir}

        SCAFTIGS={input.scaftigs}
        CUTFA={input.scaftigs_cut}
        BED={input.scaftigs_bed}
        COVERAGE={input.coverage}

        pigz -dkf $CUTFA
        pigz -dkf $COVERAGE

        set +e
        concoct \
        --threads {threads} \
        --basename {params.bin_prefix} \
        --coverage_file ${{COVERAGE%.gz}} \
        --composition_file ${{CUTFA%.gz}} \
        --clusters {params.clusters} \
        --kmer_length {params.kmer_length} \
        --length_threshold {params.length_threshold} \
        --read_length {params.read_length} \
        --total_percentage_pca {params.total_percentage_pca} \
        --seed {params.seed} \
        --iterations {params.iterations} \
        {params.no_cov_normalization} \
        {params.no_total_coverage} \
        {params.no_original_data} \
        {params.coverage_out} \
        2> {log}

        rm -rf ${{COVERAGE%.gz}}
        rm -rf ${{CUTFA%.gz}}

        cat {params.bin_prefix}_log.txt >> {log}

        grep -oEi 'Not enough contigs pass the threshold filter' {params.bin_prefix}_log.txt
        grepcode=$?
        if [ $grepcode -eq 0 ]
        then
            echo "Not enough contigs pass the threshold, touch empty concoct output directory"
        fi

        if [ -f "{params.bin_prefix}_clustering_gt{params.length_threshold}.csv" ];
        then
            merge_cutup_clustering.py \
            {params.bin_prefix}_clustering_gt{params.length_threshold}.csv \
            > {params.bin_prefix}_clustering_merged.csv

            pigz -dkf $SCAFTIGS

            extract_fasta_bins.py \
            ${{SCAFTIGS%.gz}} \
            {params.bin_prefix}_clustering_merged.csv \
            --output_path {params.mags_dir}

            rm -rf ${{SCAFTIGS%.gz}}

            python {params.wrapper_dir}/concoct_postprocess.py \
            {params.mags_dir} \
            {params.bin_prefix}

            for FILESTR in `ls {params.bin_prefix}*`
            do
                if [ -f $FILESTR ] && [ "${{FILESTR##*.}}" != "gz" ];
                then
                    pigz -f $FILESTR
                fi
            done
        fi

        touch {output}
        '''


if config["params"]["binning"]["concoct"]["do"]:
    rule binning_concoct_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "mags/{binning_group}.{assembly_group}.{assembler}/concoct/binning_done"),
                zip,
                binning_group=ASSEMBLY_GROUPS["binning_group"],
                assembly_group=ASSEMBLY_GROUPS["assembly_group"],
                assembler=ASSEMBLY_GROUPS["assembler"])

else:
    rule binning_concoct_all:
        input:


localrules:
    binning_metabat2_all,
    binning_maxbin2_all,
    binning_concoct_all