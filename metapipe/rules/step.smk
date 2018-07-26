simulate_output = expand(
    [
        "{simulate}/species_metadata.tsv", "{simulate}/merged_genome.fasta",
        "{simulate}/{out_prefix}_{read}.fastq",
        "{simulate}/{out_prefix}_abundance.txt"
    ],
    simulate=config["results"]["simulate"]["genome"],
    out_prefix=config["params"]["simulate"]["output_prefix"]["X100"],
    read=["R1", "R2"])

fastqc_output = expand(
    "{fastqc}/{sample}_{read}_fastqc.{out}",
    fastqc=config["results"]["raw"]["fastqc"],
    sample=_samples.index,
    read=["1", "2"],
    out=["html", "zip"])

trim_output = expand(
    "{trim}/{sample}.trimmed.{read}.fq.gz",
    trim=config["results"]["trim"],
    sample=_samples.index,
    read=["1", "2", "single"])

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
        "{depth}/{sample}.metabat2.depth.txt", "{logs}/{sample}.metabat2.done",
        "{logs}/{sample}.metabat2.log"
    ],
    depth=config["results"]["binning"]["depth"],
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
    [
        "{checkm}/checkm_out/{sample}.checkm.txt",
        "{checkm}/checkm_data/{sample}.checkm_out/"
    ],
    checkm=config["results"]["checkm"],
    sample=_samples.index)
'''
dereplication_output = expand(

)

classification_output = expand(

)

annotation_output = expand(

)
'''

trim_target = (fastqc_output + trim_output)
rmhost_target = (trim_target + rmhost_output)

if config["params"]["assembly"]["megahit"]["do"]:
    assembly_output = (megahit_output)
if config["params"]["assembly"]["idba_ud"]["do"]:
    assembly_output = (assembly_output + idba_ud_output)
if config["params"]["assembly"]["metaspades"]["do"]:
    assembly_output = (assembly_output + metaspades_output)

if config["params"]["rmhost"]["do"]:
    assembly_target = (rmhost_target + assembly_output)
else:
    assembly_target = (trim_target + assembly_output)

alignment_target = (assembly_target + alignment_output)

if config["params"]["binning"]["metabat2"]["do"]:
    binning_output = (metabat2_output)
if config['params']["binning"]["maxbin2"]["do"]:
    binning_output = (binning_output + maxbin2_output)

binning_target = (alignment_target + binning_output)
checkm_target = (binning_target + checkm_output)
'''
dereplication_target = (cehckm_target + dereplication_output)
classification_target = (drep_target + classification_output)
annotation_target = (classification_target + annotation_output)
'''

all_target = (fastqc_output + trim_output + rmhost_output + assembly_output +
              alignment_output + binning_output + checkm_output)
