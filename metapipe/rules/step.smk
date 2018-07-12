fastqc_output = expand(
    "{fastqc}/{sample}_{read}_fastqc.{out}",
    fastqc=config["results"]["raw"]["fastqc"],
    read=["1", "2"],
    sample=samples.keys(),
    out=["html", "zip"]
)

trim_output = expand(
    "{trim}/{sample}.trimmed.{read}.fq.gz",
    trim=config["results"]["trim"],
    sample=samples.keys(),
    read=["1", "2", "single"]
)

rmhost_output = expand(
    ["{rmhost}/{sample}.rmhost.{read}.fq.gz",
    "{rmhost}/{sample}.rmhost.flagstat"],
    rmhost=config["results"]["rmhost"],
    sample=samples.keys(),
    read=["1", "2"]
)

assembly_output = expand(
    "{assembly}/{sample}.megahit_out/{sample}.contigs.fa",
    assembly=config["results"]["assembly"],
    sample=samples.keys()
)

alignment_output = expand(
    ["{alignment}/{sample}.sorted.bam",
    "{alignment}/{sample}.flagstat"],
    alignment=config["results"]["alignment"],
    sample=samples.keys()
)

binning_output = expand(
    ["{logs}/{sample}_binning.done",
    "{logs}/{sample}.metabat2.log"],
    logs=config["logs"]["binning"],
    sample=samples.keys()
)

checkm_output = expand(
    "{checkm}/checkm_out/{sample}.checkm.txt",
    checkm=config["results"]["checkm"],
    sample=samples.keys()
)

'''
dereplication_output = expand(

)

classification_output = expand(

)

annotation_output = expand(

)
'''

'''
trim_target           = (fastqc_output + trim_output)
rmhost_target         = (trim_target + rmhost_output)
assembly_target       = (rmhost_target + assembly_output)
alignment_target      = (assembly_target + alignment_output)
binning_target        = (alignment_target + binning_output)
checkm_target         = (binning_target + checkm_output)
dereplication_target  = (cehckm_target + dereplication_output)
classification_target = (drep_target + classification_output)
annotation_target     = (classification_target + annotation_output)
'''

all_output = (
    fastqc_output +
    trim_output +
    rmhost_output +
    assembly_output +
    alignment_output +
    binning_output +
    checkm_output
)