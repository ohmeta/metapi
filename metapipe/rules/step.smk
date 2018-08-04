simulation_output = expand(
    [
        "{simulation}/species_metadata.tsv",
        "{simulation}/merged_genome.fasta",
        "{simulation}/{output_prefix}_{read}.fq.gz",
        "{simulation}/{output_prefix}_abundance.txt"
    ],
    simulation=config["results"]["simulation"]["genomes"],
    output_prefix=config["params"]["simulation"]["output_prefix"],
    read=["1", "2"])

fastqc_output = expand(
    [
        "{fastqc}/{sample}_{read}_fastqc.{out}",
        "{multiqc}/fastqc_multiqc_report.html",
        "{multiqc}/fastqc_multiqc_report_data"
    ],
    fastqc=config["results"]["raw"]["fastqc"],
    multiqc=config["results"]["raw"]["multiqc"],
    sample=_samples.index,
    read=["1", "2"],
    out=["html", "zip"])

oas1_output = expand(
    ["{trimming}/{sample}.trimmed.{read}.fq.gz",
     "{trimming}/{sample}.trimmed.stat_out"],
    trimming=config["results"]["trimming"],
    read=["1", "2", "single"],
    sample=_samples.index
)

sickle_output = expand(
    "{trimming}/{sample}.trimmed.{read}.fq.gz",
    trimming=config["results"]["trimming"],
    sample=_samples.index,
    read=["1", "2", "single"])

fastp_output = expand(
    [
        "{trimming}/{sample}.trimmed.{read}.fq.gz",
        "{trimming}/{sample}.fastp.html", "{trimming}/{sample}.fastp.json"
    ],
    sample=_samples.index,
    trimming=config["results"]["trimming"],
    read=["1", "2"])

rmhost_output = expand(
    [
        "{rmhost}/{sample}.rmhost.flagstat.txt",
        "{rmhost}/{sample}.rmhost.{read}.fq.gz"
    ],
    rmhost=config["results"]["rmhost"],
    sample=_samples.index,
    read=["1", "2"])

megahit_output = expand(
    "{assembly}/{sample}.megahit_out/{sample}.contigs.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index)

idba_ud_output = expand(
    "{assembly}/{sample}.idba_ud_out/{sample}.scaffolds.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index)

metaspades_output = expand(
    "{assembly}/{sample}.metaspades_out/{sample}.scaffolds.fa.gz",
    assembly=config["results"]["assembly"],
    sample=_samples.index)

alignment_output = expand(
    ["{alignment}/{sample}.flagstat", "{alignment}/{sample}.sorted.bam"],
    alignment=config["results"]["alignment"],
    sample=_samples.index)

metabat2_output = expand(
    [
        "{depth}/{sample}.metabat2.depth.txt", "{bins}/{sample}.metabat2_out",
        "{logs}/{sample}.metabat2.done", "{logs}/{sample}.metabat2.log"
    ],
    depth=config["results"]["binning"]["depth"],
    bins=config["results"]["binning"]["bins"],
    logs=config["logs"]["binning"]["metabat2"],
    sample=_samples.index)

maxbin2_output = expand(
    [
        "{depth}/{sample}.bbmap.depth.txt",
        "{depth}/{sample}.maxbin2.depth.txt",
        "{bins}/{sample}.maxbin2_out/{sample}.bin.summary"
    ],
    depth=config["results"]["binning"]["depth"],
    bins=config["results"]["binning"]["bins"],
    sample=_samples.index)

checkm_output = expand(
    ["{out}/{sample}.checkm.txt", "{data}/{sample}"],
    out=config["results"]["checkm"]["out"],
    data=config["results"]["checkm"]["data"],
    sample=_samples.index)
'''
dereplication_output = expand(

)

classification_output = expand(

)
'''
annotation_output = expand(
    "{prokka}/{sample}.prokka_out",
    prokka=config["results"]["annotation"]["prokka"],
    sample=_samples.index)

profilling_output = expand(
    ["{bowtie2out}/{sample}.bowtie2.gz",
     "{profile}/{sample}.metaphlan2.profile",
     "{metaphlan2}/metaphlan2.merged.profile"],
    bowtie2out=config["results"]["profilling"]["metaphlan2"]["bowtie2_out"],
    profile=config["results"]["profilling"]["metaphlan2"]["profile"],
    metaphlan2=config["results"]["profilling"]["metaphlan2"]["base_dir"],
    sample=_samples.index)

if config["params"]["trimming"]["oas1"]["do"]:
    trimming_output = (oas1_output)
if config["params"]["trimming"]["sickle"]["do"]:
    trimming_output = (sickle_output)
if config["params"]["trimming"]["fastp"]["do"]:
    trimming_output = (fastp_output)
trimming_target = (fastqc_output + trimming_output)

rmhost_target = (trimming_target + rmhost_output)

if config["params"]["assembly"]["megahit"]["do"]:
    assembly_output = (megahit_output)
if config["params"]["assembly"]["idba_ud"]["do"]:
    assembly_output = (assembly_output + idba_ud_output)
if config["params"]["assembly"]["metaspades"]["do"]:
    assembly_output = (assembly_output + metaspades_output)

if config["params"]["rmhost"]["do"]:
    assembly_target = (rmhost_target + assembly_output)
else:
    assembly_target = (trimming_target + assembly_output)

alignment_target = (assembly_target + alignment_output)

if config["params"]["binning"]["metabat2"]["do"]:
    binning_output = (metabat2_output)
if config['params']["binning"]["maxbin2"]["do"]:
    binning_output = (binning_output + maxbin2_output)

binning_target = (alignment_target + binning_output)

checkm_target = (binning_target + checkm_output)

annotation_target = (checkm_target + annotation_output)

profilling_target = (annotation_target + profilling_output)
'''
dereplication_target = (cehckm_target + dereplication_output)
classification_target = (drep_target + classification_output)
annotation_target = (classification_target + annotation_output)
'''

all_target = (
    fastqc_output + trimming_output + rmhost_output + assembly_output +
    alignment_output + binning_output + checkm_output + annotation_output)
