fastp_output = expand([
    "{trimming_dir}/trimmed.{sample}_{read}.fq.gz",
    "{trimming_dir}/trimmed.{sample}_fastp.html",
    "{trimming_dir}/{sample}_fastp.json"
],
                      trimming_dir=config["results"]["trimming"],
                      sample=_samples.index,
                      read=["1", "2"])

rmhost_output = expand([
    "{rmhost_dir}/rmhostaa.{sample}_{read}.fq.gz",
    "{rmhost_dir}/{sample}_rmhostaa_flagstat.txt",
    "{rmhost_dir}/{sample}_host_sorted.bam"
],
                       rmhost_dir=config["results"]["rmhost"],
                       sample=_samples.index,
                       read=["1", "2"])

rmhost_target = (fastp_output + rmhost_output)

