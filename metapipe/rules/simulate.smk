if config["params"]["simulate"]["do"]:
    rule genome_download:
        output:
            dir = config["results"]["simulate"]["genome"],
            metadata = os.path.join(config["results"]["simulate"]["genome"], "species_metadata.tsv")
        params:
            taxid = ",".join(config["results"]["simulate"]["taxid"])
        shell:
            "ncbi-genome-download "
            "--format fasta,assembly-report "
            "--assembly-level complete "
            "--taxid {params.taxid} "
            "--output-folder {output.dir} "
            "--human-readable "
            "--retries 3 "
            "-m {output.metadata} bacteria"

'''
rule simulate:
    input:

    output:

    shell:
        """
        iss 
        """
'''