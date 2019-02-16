fastp_output = expand([
    "{trimming}/{sample}.trimmed.{read}.fq.gz",
    "{trimming}/{sample}.fastp.html",
    "{trimming}/{sample}.fastp.json",
    "{trimming}/fastp_multiqc_report.html",
    "{trimming}/fastp_multiqc_report_data"
],
                      trimming=config["results"]["trimming"],
                      sample=_samples.id,
                      read=["1", "2"])

rmhost_output = expand([
    "{rmhost}/rmhost.{sample}.rmhost.{read}.fq.gz",
    "{rmhost}/{sample}.rmhost.flagstat.txt",
    "{rmhost}/{sample}.host.sorted.bam"
],
                       rmhost=config["results"]["rmhost"],
                       sample=_samples.index,
                       read=["1", "2"])

rmhost_target = (fastp_output + rmhost_output)

