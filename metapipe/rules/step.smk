simulate_output = expand(
    "{simulate}/species_metadata.tsv",
    simulate=config["results"]["simulate"]["genome"]
)

fastqc_output = expand(
    "{fastqc}/{sample}_{read}_fastqc.{out}",
    fastqc=config["results"]["raw"]["fastqc"],
    sample=_samples.index,
    read=["1", "2"],
    out=["html", "zip"]
)

trim_output = expand(
    "{trim}/{sample}.trimmed.{read}.fq.gz",
    trim=config["results"]["trim"],
    sample=_samples.index,
    read=["1", "2", "single"]
)

rmhost_output = expand(
    ["{rmhost}/{sample}.rmhost.flagstat.txt",
     "{rmhost}/{sample}.rmhost.{read}.fq.gz"],
    rmhost=config["results"]["rmhost"],
    sample=_samples.index,
    read=["1", "2"]
)

assembly_output = expand(
    "{assembly}/{sample}.megahit_out/{sample}.contigs.fa",
    assembly=config["results"]["assembly"],
    sample=_samples.index
)

alignment_output = expand(
    ["{alignment}/{sample}.flagstat",
     "{alignment}/{sample}.sorted.bam"],
    alignment=config["results"]["alignment"],
    sample=_samples.index
)

binning_output = expand(
    ["{binning}/coverage/{sample}.metabat2.depth.txt",
     "{logs}/{sample}.binning.done",
     "{logs}/{sample}.metabat2.log"],
    binning=config["results"]["binning"],
    logs=config["logs"]["binning"]["metabat2"],
    sample=_samples.index
)

checkm_output = expand(
    ["{checkm}/checkm_out/{sample}.checkm.txt",
     "{checkm}/checkm_data/{sample}.checkm_out/"],
    checkm=config["results"]["checkm"],
    sample=_samples.index
)

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

if config["params"]["rmhost"]["do"]:
    assembly_target = (rmhost_target + assembly_output)
else:
    assembly_target = (trim_target + assembly_output)

alignment_target = (assembly_target + alignment_output)
binning_target = (alignment_target + binning_output)
checkm_target = (binning_target + checkm_output)

'''
dereplication_target = (cehckm_target + dereplication_output)
classification_target = (drep_target + classification_output)
annotation_target = (classification_target + annotation_output)
'''

all_target = (
    fastqc_output +
    trim_output +
    rmhost_output +
    assembly_output +
    alignment_output +
    binning_output +
    checkm_output
)
