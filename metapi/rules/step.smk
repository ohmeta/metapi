sra2fq_output = expand(
    "{sra2fq}/{sample}{read}.fq.gz",
    sra2fq=config["results"]["sra2fq"],
    read=[".1", ".2"] if IS_PE else "",
    sample=_samples.index.unique())

simulation_output = expand([
    "{simulation}/species_metadata.tsv", "{simulation}/{sample}_genome.fa",
    "{raw}/{sample}{read}.fq.gz", "{raw}/{sample}_abundance.txt"
],
                           simulation=config["results"]["simulation"],
                           raw=config["results"]["raw"]["reads"],
                           read=[".1", ".2"] if IS_PE else "",
                           sample=_samples.index.unique())

fastqc_output = expand([
    "{fastqc}/{sample}/done",
    "{multiqc}/fastqc_multiqc_report.html",
    "{multiqc}/fastqc_multiqc_report_data"
],
                       fastqc=config["results"]["raw"]["fastqc"],
                       multiqc=config["results"]["raw"]["multiqc"],
                       sample=_samples.index.unique())

oas1_output = expand([
    "{trimming}/{sample}.trimmed{read}.fq.gz",
    "{trimming}/{sample}.trimmed.stat_out"
],
                     trimming=config["results"]["trimming"],
                     read=[".1", ".2", ".single"] if IS_PE else "",
                     sample=_samples.index.unique())

sickle_output = expand(
    "{trimming}/{sample}.trimmed{read}.fq.gz",
    trimming=config["results"]["trimming"],
    sample=_samples.index.unique(),
    read=[".1", ".2", ".single"] if IS_PE else "")

fastp_output = expand([
#    temp("{trimming}/{sample}.trimmed{read}.fq.gz"),
    "{trimming}/{sample}.fastp.html",
    "{trimming}/{sample}.fastp.json",
    "{trimming}/fastp_multiqc_report.html",
    "{trimming}/fastp_multiqc_report_data"
],
                      sample=_samples.index.unique(),
                      trimming=config["results"]["trimming"],
                      read=[".1", ".2"] if IS_PE else "")

rmhost_output = expand([
    "{rmhost}/{sample}.rmhost.flagstat.txt",
    "{rmhost}/{sample}.rmhost{read}.fq.gz"
],
                       rmhost=config["results"]["rmhost"],
                       sample=_samples.index.unique(),
                       read=[".1", ".2"] if IS_PE else "")

raw_report_output = expand(
    "{reportout}/raw.stats.tsv",
    reportout=config["results"]["report"]["base_dir"])

trimming_report_output = expand(
    "{reportout}/trimming.stats.tsv",
    reportout=config["results"]["report"]["base_dir"])

rmhost_report_output = expand(
    "{reportout}/rmhost.stats.tsv",
    reportout=config["results"]["report"]["base_dir"])

qc_report_output = expand(
    "{reportout}/qc.stats.tsv",
    reportout=config["results"]["report"]["base_dir"])

megahit_output = expand(
    "{assembly}/{sample}.megahit_out/{sample}.megahit.scaftigs.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index.unique())

idba_ud_output = expand(
    "{assembly}/{sample}.idba_ud_out/{sample}.idba_ud.scaftigs.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index.unique())

metaspades_output = expand(
    "{assembly}/{sample}.metaspades_out/{sample}.metaspades.scaftigs.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index.unique())

spades_output = expand(
    "{assembly}/{sample}.spades_out/{sample}.spades.scaftigs.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index.unique())

assembly_report_output = expand(
    [
        "{assembly}/{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.seqtk.comp.tsv.gz",
        "{reportout}/{assembler}.assembly.summary.tsv"
    ],
    assembly=config["results"]["assembly"],
    sample=_samples.index.unique(),
    assembler=config["params"]["assembler"],
    reportout=config["results"]["report"]["base_dir"]
)

upload_cnsa_output = expand(
    "{upload}/{target}.xlsx",
    upload=config["results"]["upload"],
    target=["MIxS_Samples", "Experiment_Run", "Genome_Assembly"]
)

coassembly_megahit_output = expand(
    "{coassembly_megahit}/final.contigs.fa.gz",
    coassembly_megahit=config["results"]["coassembly"]["megahit"])

demultiplex_kraken2_output = expand(
    [
        "{demultiplex_kraken2}/{sample}.demultiplex_out/{sample}.demultiplex.done",
        "{demultiplex_kraken2}/merged"
    ],
    demultiplex_kraken2=config["results"]["coassembly"]["demultiplex_kraken2"],
    sample = _samples.index.unique())

metaquast_output = expand([
    "{metaquast}/{sample}.{assembler}.metaquast_out/report.html",
    "{metaquast}/{sample}.{assembler}.metaquast_out/icarus.html",
    "{metaquast}/{sample}.{assembler}.metaquast_out/combined_reference/report.tsv",
    "{metaquast}/{sample}.{assembler}.metaquast_out/icarus_viewers",
    "{metaquast}/{sample}.{assembler}.metaquast_out/krona_charts",
    "{metaquast}/{sample}.{assembler}.metaquast_out/not_aligned",
    "{metaquast}/{sample}.{assembler}.metaquast_out/runs_per_reference",
    "{metaquast}/{sample}.{assembler}.metaquast_out/summary",
    "{metaquast}/metaquast_multiqc_report.html",
    "{metaquast}/metaquast_multiqc_report_data"
],
                          metaquast=config["results"]["metaquast"],
                          assembler=config["params"]["assembler"],
                          sample=_samples.index.unique())

prodigal_output = expand([
    "{prediction}/{sample}.prodigal_out/{sample}.{assembler}.pep.faa",
    "{prediction}/{sample}.prodigal_out/{sample}.{assembler}.cds.ffn",
    "{prediction}/{sample}.prodigal_out/{sample}.{assembler}.cds.gff",
    "{prediction}/{sample}.prodigal_out/{sample}.{assembler}.score.gff"
],
                          prediction=config["results"]["prediction"],
                          assembler=config["params"]["assembler"],
                          sample=_samples.index.unique())

alignment_bwa_output = expand([
    "{alignment}/{sample}.bwa_out/{sample}.{assembler}.flagstat",
#    "{alignment}/{sample}.bwa_out/{sample}.{assembler}.sorted.bam",
#    "{alignment}/{sample}.bwa_out/{sample}.{assembler}.sorted.bam.bai",
    "{alignment}/scaftigs_flagstat_summary.tsv"
],
                          alignment=config["results"]["alignment"],
                          assembler=config["params"]["assembler"],
                          sample=_samples.index.unique())

metabat2_output = expand([
    "{depth}/{sample}.{assembler}.metabat2.depth.txt",
    "{bins}/{sample}.{assembler}.metabat2_out",
    "{logs}/{sample}.{assembler}.metabat2.log"
],
                         depth=config["results"]["binning"]["depth"],
                         bins=config["results"]["binning"]["bins"],
                         logs=config["logs"]["binning"]["metabat2"],
                         assembler=config["params"]["assembler"],
                         sample=_samples.index.unique())

maxbin2_output = expand([
    "{depth}/{sample}.{assembler}.bbmap.depth.txt",
    "{depth}/{sample}.{assembler}.maxbin2.depth.txt",
    "{bins}/{sample}.{assembler}.maxbin2_out/{sample}.bin.summary"
],
                        depth=config["results"]["binning"]["depth"],
                        bins=config["results"]["binning"]["bins"],
                        assembler=config["params"]["assembler"],
                        sample=_samples.index.unique())

cobin_prediction_output = expand([
    "{cds}/{sample}/{sample}.{assembler}.cds.fa.gz",
    "{cds}/{sample}/{sample}.{assembler}.cds.gff.gz"
],
    cds=config["results"]["cobinning"]["cds"],
    sample=_samples.index.unique(),
    assembler=config["params"]["assembler"])

cobin_vsearch_clust_output = expand([
    "{cds}/{sample}/{sample}.{assembler}.cds.uc.gz",
    "{cds}/{sample}/{sample}.{assembler}.cds.marker.fa.gz"
],
    cds=config["results"]["cobinning"]["cds"],
    sample=_samples.index.unique(),
    assembler=config["params"]["assembler"])

cobin_alignment_cds_output = expand(
    "{depth}/{sample_}/{sample_}.{sample}.{assembler}.{binner}.depth.txt.gz",
    depth=config["results"]["cobinning"]["depth"],
    sample_= _samples_id,
    sample=_samples.index.unique(),
    assembler=config["params"]["assembler"],
    binner=config["params"]["binning"]["binner"])

checkm_lineage_wf_output = expand([
    "{out}/{sample}.{assembler}.{binner}.checkm.txt",
    "{data}/{sample}.{assembler}.{binner}.checkm.data.tar.gz"
],
                        out=config["results"]["checkm"]["out"],
                        data=config["results"]["checkm"]["data"],
                        assembler=config["params"]["assembler"],
                        binner=config["params"]["binning"]["binner"],
                        sample=_samples.index.unique())

checkm_report_output = expand(
    "{basedir}/{assembler}.{binner}.checkm.out.tsv",
    basedir=config["results"]["checkm"]["base_dir"],
    assembler=config["params"]["assembler"],
    binner=config["params"]["binning"]["binner"]
)

checkm_link_bins_output = expand(
    [
        "{basedir}/bins.{assembler}.{binner}_out.hq",
        "{basedir}/bins.{assembler}.{binner}_out.mq",
        "{basedir}/bins.{assembler}.{binner}_out.lq",
        "{basedir}/bins.{assembler}.{binner}_out.hmq"
    ],
    basedir=config["results"]["checkm"]["base_dir"],
    assembler=config["params"]["assembler"],
    binner=config["params"]["binning"]["binner"]
)

'''
checkm_coverage_output = expand(
    "{coverage}/{sample}.{assembler}.checkm_coverage.tsv",
    coverage=config["results"]["checkm"]["coverage"],
    assembler=config["params"]["assembler"],
    sample=_samples.index)

checkm_profile_output = expand(
    "{profile}/{sample}.{assembler}.checkm_profile.tsv",
    profile=config["results"]["checkm"]["profile"],
    assembler=config["params"]["assembler"],
    sample=_samples.index)

dereplication_output = expand()
'''

prokka_output = expand(
    [
        "{prokka}/{sample}.{assembler}.prokka_out/done",
        "{multiqc_prokka}/prokka_multiqc_report.html",
        "{multiqc_prokka}/prokka_multiqc_report_data"
    ],
    prokka=config["results"]["annotation"]["prokka"],
    multiqc_prokka=config["results"]["annotation"]["multiqc_prokka"],
    assembler=config["params"]["assembler"],
    sample=_samples.index.unique())

kraken2_output = expand(
    [
        "{kraken2}/{sample}.kraken2.output",
        "{kraken2}/{sample}.kraken2.report"
    ],
    kraken2=config["results"]["classification"]["kraken2"],
    sample=_samples.index.unique())

gtdbtk_output = expand(
    "{gtdbtkout}/hmq.bins.{assembler}.{binner}.gtdbtk_out",
    gtdbtkout=config["results"]["classification"]["gtdbtk"]["base_dir"],
    assembler=config["params"]["assembler"],
    binner=config["params"]["binning"]["binner"]
)

metaphlan2_profiling_output = expand(
    [
        "{bowtie2out}/{sample}.bowtie2.bz2",
        "{profile}/{sample}.metaphlan2.profile",
        "{metaphlan2}/metaphlan2.merged.profile"
    ],
    bowtie2out=config["results"]["profiling"]["metaphlan2"]["bowtie2_out"],
    profile=config["results"]["profiling"]["metaphlan2"]["profile"],
    metaphlan2=config["results"]["profiling"]["metaphlan2"]["base_dir"],
    sample=_samples.index.unique())

jgi_profiling_output = expand(
    "{depth}/{sample}.jgi.depth.gz",
    depth=config["results"]["profiling"]["jgi"]["depth"],
    sample=_samples.index.unique())

jgi_profile_merge_output = expand(
    [
        "{profile}/abundance_profile.tsv",
        "{profile}/depth_profile.tsv",
        "{profile}/abundance_profile_{level}.tsv",
    ],
    profile=config["results"]["profiling"]["jgi"]["profile"],
    level=["superkingdom", "phylum", "order", "class", "family", "genus", "species", "strain"]
)

humann2_profiling_output = expand(
    "{humann2}/{sample}.humann2_out/{sample}_{target}.tsv",
    humann2=config["results"]["profiling"]["humann2"],
    sample=_samples.index.unique(),
    target=["genefamilies", "pathabundance", "pathcoverage"])

humann2_postprocess_output = expand(
    [
        "{humann2}/{sample}.humann2_out/{sample}_{target}.{norm}.tsv",
        "{humann2}/{sample}.humann2_out/{sample}_group-{group}-profile.tsv"
    ],
    humann2=config["results"]["profiling"]["humann2"],
    sample=_samples.index.unique(),
    target=["genefamilies", "pathabundance", "pathcoverage"],
    norm=config["params"]["profiling"]["humann2"]["normalize_method"],
    group=config["params"]["profiling"]["humann2"]["map_database"])

humann2_join_split_output = expand(
    [
        "{humann2}/{target}_joined.tsv",
        "{humann2}/{target}_joined_{suffix}.tsv"
    ],
    humann2=config["results"]["profiling"]["humann2"],
    target=["gene_families", "path_abundance", "path_coverage"] + \
           config["params"]["profiling"]["humann2"]["map_database"],
    suffix=["straified", "unstraified"])

burst_output = expand(
    "{burst}/{sample}.reads.burst.b6",
    burst=config["results"]["burst"],
    sample=_samples.index.unique())

#----------------------------------------------------------------#

qc_output = ([])
if config["params"]["fastqc"]["do"]:
    qc_output = (fastqc_output)
qc_target = (qc_output)

trimming_output = ([])
if config["params"]["trimming"]["oas1"]["do"]:
    trimming_output = (oas1_output)
if config["params"]["trimming"]["sickle"]["do"]:
    trimming_output = (sickle_output)
if config["params"]["trimming"]["fastp"]["do"]:
    trimming_output = (fastp_output)
trimming_target = (qc_target + trimming_output)

rmhost_output_ = ([])
if config["params"]["rmhost"]["do"]:
    rmhost_output_ = (rmhost_output)
rmhost_target = (trimming_target + rmhost_output_)

assembly_output = ([])
if config["params"]["assembly"]["megahit"]["do"]:
    assembly_output = (megahit_output)
if config["params"]["assembly"]["idba_ud"]["do"]:
    assembly_output = (assembly_output + idba_ud_output)
if config["params"]["assembly"]["metaspades"]["do"]:
    assembly_output = (assembly_output + metaspades_output)
if config["params"]["assembly"]["spades"]["do"]:
    assembly_output = (assembly_output + spades_output)

if config["params"]["coassembly"]["megahit"]["do"]:
    assembly_output = (assembly_output + coassembly_megahit_output)

if config["params"]["coassembly"]["demultiplex_kraken2"]["do"]:
    assembly_output = (assembly_output + demultiplex_kraken2_output)

assembly_target = ([])
if config["params"]["begin"] == "assembly":
    if config["params"]["reads_format"] == "fastq":
        assembly_target = (assembly_output)
    elif config["params"]["reads_format"] == "sra":
        assembly_target = (sra2fq_output + assembly_output)
elif config["params"]["rmhost"]["do"]:
    assembly_target = (rmhost_target + assembly_output)
    if config["params"]["report"]["seqkit"]["do"]:
        assembly_target = (assembly_target +
                           raw_report_output +
                           trimming_report_output +
                           rmhost_report_output +
                           qc_report_output)
else:
    assembly_target = (trimming_target + assembly_output)
    if config["params"]["report"]["seqkit"]["do"]:
        assembly_target = (assembly_target +
                           raw_report_output +
                           trimming_report_output)

if config["params"]["assembly"]["report"]["do"]:
    assembly_target = (assembly_target + assembly_report_output)

if config["params"]["metaquast"]["do"]:
    assembly_target = (assembly_target + metaquast_output)

prediction_output = ([])
if config["params"]["prediction"]["prodigal"]["do"]:
    prediction_output = (prodigal_output)
prediction_target = (assembly_target + prediction_output)

upload_output = ([])
if config["upload"]["do"]:
    upload_output = (upload_cnsa_output)
upload_target = (prediction_target + upload_output)

alignment_output = ([])
if config["params"]["alignment"]["do"]:
    alignment_output = (alignment_bwa_output)
alignment_target = (upload_target + alignment_output)

binning_output = ([])
if config["params"]["binning"]["metabat2"]["do"]:
    binning_output = (metabat2_output)
if config["params"]["binning"]["maxbin2"]["do"]:
    binning_output = (binning_output + maxbin2_output)
binning_target = (alignment_target + binning_output)

cobinning_output = ([])
if config["params"]["cobinning"]["do"]:
    cobinning_output = (cobin_prediction_output +
                        cobin_vsearch_clust_output +
                        cobin_alignment_cds_output)
binning_target = (binning_target + cobinning_output)

checkm_output = ([])
if config["params"]["checkm"]["do"]:
    checkm_output = (checkm_lineage_wf_output)
    # checkm_output = (checkm_lineage_wf_output +
    #                  checkm_coverage_output +
    #                  checkm_profile_output)
if config["params"]["checkm"]["report"]:
    checkm_output = (checkm_output +
                     checkm_report_output +
                     checkm_link_bins_output)
checkm_target = (binning_target + checkm_output)

annotation_output = ([])
if config["params"]["annotation"]["prokka"]["do"]:
    annotation_output = (prokka_output)
annotation_target = (checkm_target + annotation_output)

classification_output = ([])
if config["params"]["classification"]["kraken2"]["do"]:
    classification_output = (kraken2_output)
if config["params"]["classification"]["gtdbtk"]["do"]:
    classification_output = (classification_output + gtdbtk_output)
classification_target = (annotation_target + classification_output)

profiling_output = ([])
if config["params"]["profiling"]["metaphlan2"]["do"]:
    profiling_output = (metaphlan2_profiling_output)
if config["params"]["profiling"]["jgi"]["do"]:
    profiling_output = (profiling_output +
                        jgi_profiling_output +
                        jgi_profile_merge_output)
if config["params"]["profiling"]["humann2"]["do"]:
    profiling_output = (profiling_output +
                        humann2_profiling_output +
                        humann2_postprocess_output +
                        humann2_join_split_output)
profiling_target = (classification_target + profiling_output)

burst_output_ = ([])
if config["params"]["burst"]["do"]:
    burst_output_ = (burst_output)
burst_target = (profiling_target + burst_output_)

all_target = burst_target

debug_target = (
    qc_output +
    trimming_output +
    rmhost_output +
    assembly_output +
    prediction_output +
    upload_output +
    alignment_output +
    binning_output +
    checkm_output +
    annotation_output +
    classification_output +
    profiling_output +
    burst_output)

# all_target == debug_target
