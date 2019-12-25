rule rmhost_md5:
    input:
        expand(os.path.join(config["results"]["rmhost"], "{{sample}}.rmhost{read}.fq.gz"),
               read=[".1", ".2"] if IS_PE else "")
    output:
        os.path.join(config["results"]["rmhost"], "{{sample}}.rmhost.fq.gz.md5")
    params:
        rmhost_dir = os.path.join(config["results"]["rmhost"], "")
    shell:
        '''
        md5sum {input} | sed 's#{params.rmhost_dir}##g' > {output}
        '''


rule assembly_md5:
    input:
        os.path.join(config["results"]["assembly"],
                     "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz"),
    output:
        os.path.join(config["results"]["assembly"],
                     "{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz.md5")
    params:
        assembly_dir = os.path.join(config["results"]["assembly"], "{sample}.{assembler}_out/")
    shell:
        '''
        md5sum {intput} | sed 's#{params.assembly_dir}##g' > {output}
        '''

       
rule generate_samples_info:
    input:
        config["params"]["samples"]
    output:
        os.path.join(config["results"]["upload"], "MIxS_Samples.xlsx")
    run:
        uploader.gen_samples_info(_samples, output, config)
       

rule generate_run_info:
    input:
        expand("{rmhost}/{sample}.rmhost{read}.fq.gz.md5",
               rmhost=config["results"]["rmhost"],
               sample=_samples.index.unique(),
               read = [".1", ".2"] if IS_PE else "")
    output:
        os.path.join(config["results"]["upload"], "Experiment_Run.xlsx")
    threads:
        config["upload"]["threads"]
    run:
        uploader.gen_info(input, output, config, threads, "sequencing_run")


rule generate_assembly_info:
    input:
        expand("{assembly}/{sample}.{assembler}_out/{sample}.{assembler}.scaftigs.fa.gz.md5",
               assembly=config["results"]["assembly"],
               sample=_samples.index.unique(),
               assembler=config["params"]["assembler"])
    output:
        os.path.join(config["results"]["upload"], "Genome_Assembly.xlsx")
    threads:
        config["upload"]["threads"]
    run:
        uploader.gen_info(input, output, config, threads, "assembly")
