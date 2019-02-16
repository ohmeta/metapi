fastp_output = expand([
    "{trimming}/{sample}_trimmed_{read,[12]}.fq.gz",
    "{trimming}/{sample}_fastp.html",
    "{trimming}/{sample}_fastp.json",
    "{trimming}/fastp_multiqc_report.html",
    "{trimming}/fastp_multiqc_report_data"
],
                      trimming=config["results"]["trimming"],
                      sample=_samples.id,
                      read=["1", "2"])

rmhost_output = expand([
    "{rmhost}/{sample}_rmhost_{read,[12]}.fq.gz",
    "{rmhost}/{sample}_rmhost_flagstat.txt",
    "{rmhost}/{sample}_host_sorted.bam"
],
                       rmhost=config["results"]["rmhost"],
                       sample=_samples.index,
                       read=["1", "2"])

rmhost_target = (fastp_output + rmhost_output)

