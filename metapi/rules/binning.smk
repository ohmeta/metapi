if config["params"]["binning"]["metabat2"]["do"]:
    rule binning_metabat2_coverage:
        input:
            bam = lambda wildcards: expand(os.path.join(
                config["output"]["alignment"],
                "bam/{{assembly_group}}.{{assembler}}.out/{sample}.align2scaftigs.sorted.bam"),
                sample=metapi.get_samples_id_by_assembly_group(SAMPLES, wildcards.assembly_group)),
            bai = lambda wildcards: expand(os.path.join(
                config["output"]["alignment"],
                "bam/{{assembly_group}}.{{assembler}}.out/{sample}.align2scaftigs.sorted.bam.bai"),
                sample=metapi.get_samples_id_by_assembly_group(SAMPLES, wildcards.assembly_group))
        output:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.coverage")
        priority:
            28
        conda:
            config["envs"]["metabat2"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{assembly_group}.{assembler}.metabat2.coverage.log")
        params:
            percent_identity = config["params"]["binning"]["metabat2"]["percent_identity"],
            min_map_qual = config["params"]["binning"]["metabat2"]["min_map_qual"],
            output_paired_contigs = "--pairedContigs %s" % \
                os.path.join(
                    config["output"]["binning"],
                    "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.paired_contigs") \
                    if config["params"]["binning"]["metabat2"]["output_paired_contigs"] \
                       else "",
            output_gc = "--outputGC %s" % \
                os.path.join(
                    config["output"]["binning"],
                    "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.gc") \
                    if config["params"]["binning"]["metabat2"]["output_gc"] \
                       else "",
            output_gc_window = "--gcWindow %s" % \
                os.path.join(
                    config["output"]["binning"],
                    "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.gc_window") \
                    if config["params"]["binning"]["metabat2"]["output_gc_window"] \
                       else "",
            output_dir = os.path.join(config["output"]["binning"],
                                      "coverage/{assembly_group}.{assembler}.out")
        shell:
            '''
            jgi_summarize_bam_contig_depths \
            --outputDepth {output.coverage} \
            --percentIdentity {params.percent_identity} \
            --minMapQual {params.min_map_qual} \
            {params.output_paired_contigs} \
            {params.output_gc} \
            {params.output_gc_window} \
            {input.bam} \
            2> {log}
            '''


    rule binning_metabat2_coverage_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.coverage"),
                   assembler=ASSEMBLERS,
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)


    rule binning_metabat2:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz"),
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.coverage")
        output:
            os.path.join(config["output"]["binning"],
                         "bins/{assembly_group}.{assembler}.out/metabat2/binning_done")
        priority:
            30
        conda:
            config["envs"]["metabat2"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{assembly_group}.{assembler}.metabat2.binning.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/metabat2/{assembly_group}.{assembler}.metabat2.benchmark.txt")
        params:
            bins_dir = os.path.join(config["output"]["binning"], "bins/{assembly_group}.{assembler}.out/metabat2"),
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/metabat2/{assembly_group}.{assembler}.metabat2.bin"),
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
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            rm -rf {params.bins_dir}

            metabat2 \
            --inFile {input.scaftigs} \
            --abdFile {input.coverage} \
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
            --verbose > {log}

            touch {output}
            '''


    rule binning_metabat2_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}.{assembler}.out/metabat2/binning_done"),
                assembler=ASSEMBLERS,
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

            #rules.binning_metabat2_coverage_all.input,
            #rules.alignment_all.input,
            #rules.assembly_all.input
       
else:
    rule binning_metabat2_all:
        input:



if config["params"]["binning"]["maxbin2"]["do"]:
    rule binning_maxbin2_coverage:
        input:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.metabat2.coverage")
        output:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.maxbin2.coverage")
        priority:
            30
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{assembly_group}.{assembler}.maxbin2.coverage.log")
        shell:
            '''
            cut -f1,3 {input.coverage} | tail -n +2 > {output.coverage}
            '''


    rule binning_maxbin2:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz"),
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.maxbin2.coverage")
        output:
            os.path.join(config["output"]["binning"],
                         "bins/{assembly_group}.{assembler}.out/maxbin2/binning_done")
        priority:
            30
        conda:
            config["envs"]["maxbin2"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{assembly_group}.{assembler}.maxbin2.binning.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/maxbin2/{assembly_group}.{assembler}.maxbin2.benchmark.txt")
        params:
            wrapper_dir = WRAPPER_DIR,
            bins_dir = os.path.join(config["output"]["binning"],
                                    "bins/{assembly_group}.{assembler}.out/maxbin2"),
            bin_prefix = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/maxbin2/{assembly_group}.{assembler}.maxbin2.bin"),
            min_contig = config["params"]["binning"]["maxbin2"]["min_contig"],
            max_iteration = config["params"]["binning"]["maxbin2"]["max_iteration"],
            prob_threshold = config["params"]["binning"]["maxbin2"]["prob_threshold"],
            plotmarker = "-plotmarker" if config["params"]["binning"]["maxbin2"]["plotmarker"] \
                else "",
            markerset = config["params"]["binning"]["maxbin2"]["markerset"]
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            rm -rf {params.bins_dir}
            mkdir -p {params.bins_dir}

            set +e
            run_MaxBin.pl \
            -thread {threads} \
            -contig {input.scaftigs} \
            -abund {input.coverage} \
            -min_contig_length {params.min_contig} \
            -max_iteration {params.max_iteration} \
            -prob_threshold {params.prob_threshold} \
            {params.plotmarker} \
            -markerset {params.markerset} \
            -out {params.bin_prefix} \
            > {log} 2>&1

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
            {params.bins_dir}

            touch {output}
            '''


    rule binning_maxbin2_all:
        input:
            expand(
                os.path.join(
                    config["output"]["binning"],
                    "bins/{assembly_group}.{assembler}.out/maxbin2/binning_done"),
                assembler=ASSEMBLERS,
                assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

            #rules.alignment_all.input,
            #rules.assembly_all.input

else:
    rule binning_maxbin2_all:
        input:


'''
if config["params"]["binning"]["canopy"]["do"]:
    rule binning_canopy_coverage:
        input:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.metabat2.coverage")
        output:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{sample}.{assembler}.out/{sample}.{assembler}.canopy.coverage")
        priority:
            30
        run:
            import pandas as pd

            df = pd.read_csv(input.coverage, sep='\t')
            df.iloc[:, [0, 3]].to_csv(output.coverage, header=None, sep='\t', index=False)
'''


if config["params"]["binning"]["concoct"]["do"]:
    rule binning_concoct_cut_bed:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa.gz"),
        output:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa"),
            scaftigs_cut = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.cut.fa"),
            scaftigs_bed = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.cut.bed")
        threads:
            1
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{assembly_group}.{assembler}.concoct.cut_bed.log")
        conda:
            config["envs"]["concoct"]
        params:
            chunk_size = config["params"]["binning"]["concoct"]["chunk_size"],
            overlap_size = config["params"]["binning"]["concoct"]["overlap_size"]
        shell:
            '''
            pigz -p {threads} -k -d -c {input.scaftigs} > {output.scaftigs}

            cut_up_fasta.py \
            {output.scaftigs} \
            --chunk_size {params.chunk_size} \
            --overlap_size {params.overlap_size} \
            --merge_last \
            --bedfile {output.scaftigs_bed} \
            > {output.scaftigs_cut} 2>{log}
            '''


    rule binning_concoct_cut_bed_all:
        input:
            expand(os.path.join(config["output"]["assembly"],
                   "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.{results}"),
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST,
                   assembler=ASSEMBLERS,
                   results=["scaftigs.fa", "scaftigs.cut.fa", "scaftigs.cut.bed"])


    rule binning_concoct_coverage:
        input:
            scaftigs_bed = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.cut.bed"),
            bam = lambda wildcards: expand(os.path.join(
                config["output"]["alignment"],
                "bam/{{assembly_group}}.{{assembler}}.out/{sample}.align2scaftigs.sorted.bam"),
                sample=metapi.get_samples_id_by_assembly_group(SAMPLES, wildcards.assembly_group)),
            bai = lambda wildcards: expand(os.path.join(
                config["output"]["alignment"],
                "bam/{{assembly_group}}.{{assembler}}.out/{sample}.align2scaftigs.sorted.bam.bai"),
                sample=metapi.get_samples_id_by_assembly_group(SAMPLES, wildcards.assembly_group))
        output:
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.concoct.coverage")
        priority:
            30
        conda:
            config["envs"]["concoct"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/coverage/{assembly_group}.{assembler}.concoct.coverage.log")
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            concoct_coverage_table.py \
            {input.scaftigs_bed} \
            {input.bam} \
            > {output.coverage} 2> {log}
            '''


    rule binning_concoct:
        input:
            scaftigs = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.fa"),
            scaftigs_cut = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.cut.fa"),
            scaftigs_bed = os.path.join(
                config["output"]["assembly"],
                "scaftigs/{assembly_group}.{assembler}.out/{assembly_group}.{assembler}.scaftigs.cut.bed"),
            coverage = os.path.join(
                config["output"]["binning"],
                "coverage/{assembly_group}.{assembler}.out/{assembly_group}.concoct.coverage")
        output:
            os.path.join(config["output"]["binning"],
                         "bins/{assembly_group}.{assembler}.out/concoct/binning_done")
        conda:
            config["envs"]["concoct"]
        log:
            os.path.join(config["output"]["binning"],
                         "logs/binning/{assembly_group}.{assembler}.concoct.binning.log")
        benchmark:
            os.path.join(config["output"]["binning"],
                         "benchmark/concoct/{assembly_group}.{assembler}.concoct.benchmark.txt")
        params:
            wrapper_dir = WRAPPER_DIR,
            bins_dir = os.path.join(config["output"]["binning"],
                                    "bins/{assembly_group}.{assembler}.out/concoct"),
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
            basename = os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/concoct/{assembly_group}.{assembler}.concoct.bin"),
        threads:
            config["params"]["binning"]["threads"]
        shell:
            '''
            rm -rf {params.bins_dir}
            mkdir -p {params.bins_dir}

            set +e

            concoct \
            --threads {threads} \
            --basename {params.basename} \
            --coverage_file {input.coverage} \
            --composition_file {input.scaftigs_cut} \
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

            cat {params.basename}_log.txt >> {log}

            grep -oEi 'Not enough contigs pass the threshold filter' {params.basename}_log.txt
            grepcode=$?
            if [ $grepcode -eq 0 ]
            then
                echo "Not enough contigs pass the threshold, touch empty concoct output directory"
            fi

            if [ -f "{params.basename}_clustering_gt{params.length_threshold}.csv" ];
            then
                merge_cutup_clustering.py \
                {params.basename}_clustering_gt{params.length_threshold}.csv \
                > {params.basename}_clustering_merged.csv

                extract_fasta_bins.py \
                {input.scaftigs} \
                {params.basename}_clustering_merged.csv \
                --output_path {params.bins_dir}

                python {params.wrapper_dir}/concoct_postprocess.py \
                {params.bins_dir} \
                {params.basename}
            fi

            touch {output}
            '''


    rule binning_concoct_all:
        input:
            expand(os.path.join(
                config["output"]["binning"],
                "bins/{assembly_group}.{assembler}.out/concoct/binning_done"),
                   assembler=ASSEMBLERS,
                   assembly_group=SAMPLES_ASSEMBLY_GROUP_LIST)

            #rules.alignment_all.input,
            #rules.assembly_all.input

else:
    rule binning_concoct_all:
        input:


localrules:
    binning_metabat2_all,
    binning_metabat2_coverage_all,
    binning_maxbin2_all,
    binning_concoct_cut_bed_all,
    binning_concoct_all