checkpoint get_bins:
    input:
        fasta = "/vamb/xxx.megahit.combined.scaftigs.fa.gz",
        clusters = "/vamb/xxx.megahit.vamb.out/clusters.tsv"

    output:
        fna = "resultdir/vamb_bins.{i}.fna"

    params:
        phamb_utils = "mag_annotation/scripts/get_contigs_vamb_bins.py", # revised script
        resultdir = lambda wildcards, output: Path(output.fna).parent,
        min_bin_size = 5000

	shell:
	    '''
	    python {params.phamb_utils} {input.fasta} {input.clusters} {params.resultdir} -m {params.min_bin_size}
	    '''

# output: resultdir/vamb_bins_all_ge5000_without_running/vamb_bins.*.fna
# each portion contains 100000 bins at max


bin_files = os.listdir("resultdir/vamb_bins_all_ge5000_without_running/.*fna")

PARTS = list(range(1,11))
PARTS = PARTS[:len(bin_files)]

# with reference to https://evodify.com/snakemake-checkpoint-tutorial/
# and https://github.com/snakemake/snakemake/issues/1122

def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the get_bins step
    '''
    checkpoint_output = checkpoints.get_bins.get(**wildcards).output[0]
    return expand('resultdir/vamb_bins.{i}.fna',
           i=PARTS)


rule dbscan_swa:
    input:
        # bins_fasta = expand("resultdir/vamb_bins_all_ge5000_without_running/vamb_bins.{part}.fna", part=PARTS)
        bins_fasta = aggregate_input
    output:
        fna = "output.fna"
        faa = "output.faa"

    params:
        dbscan_swa_script = "/home/yepeng/software/DBSCAN-SWA-1/bin/dbscan-swa.py", # revised script
        prefix = "test",
        output_dir = expand("dbscan_swa/vamb_bins.{part}", part=PARTS)

    threads:
	    20
	    # config["params"]["dbscan_swa"]["threads"]

    shell:
        '''
        # for loop here?
        
        python {params.dbscan_swa_script} --input {input.bins_fasta} --output {params.output_dir} --thread_num {threads} --prefix {params.prefix} --add_annotation none
        cat {params.output_dir}/*.fna > {output.fna}
        cat {params.output_dir}/*.faa > {output.faa}
        head -1 {params.output_dir}/*1*summary.txt > {sample}_DBSCAN_SWA_summary_all.txt
        tail -n +2 -q {params.output_dir}/*prophage_summary.txt >> {sample}_DBSCAN_SWA_summary_all.txt
        '''
# result file: output_dir/*summary.txt output_dir/*.fna output_dir/*.faa
