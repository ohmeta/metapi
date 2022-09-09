rule dereplicate_gene_prepare:
    input:
        expand(os.path.join(
            config["output"]["predict"],
            "scaftigs_gene/{binning_group}.{assembly_group}.{{assembler}}.prodigal/{binning_group}.{assembly_group}.{{assembler}}.prodigal.ffn"),
            zip,
            binning_group=ASSEMBLY_GROUP["binning_group"],
            assembly_group=ASSEMBLY_GROUP["assembly_group"])
    output:
        ffn = os.path.join(config["output"]["predict"],
                           "scaftigs_gene_merged/{assembler}.prodigal.scaftigs.gene.merged.ffn"),
        metadata = os.path.join(config["output"]["predict"],
                                "scaftigs_gene_merged/{assembler}.prodigal.scaftigs.gene.merged.ffn.metadata")
    run:
        from Bio import SeqIO

        mg_count = 0

        with open(output.ffn, 'w') as fh, open(output.metadata, 'w') as mh:
            mh.write("mg_id\tcds_id\tmg_name\tcds_name\n")
            for i in input:
                mg_count += 1
                cds_count = 0
                for seq_record in SeqIO.parse(i, "fasta"):
                    cds_count += 1
                    mh.write(
                        f"MG_{mg_count}\tCDS_{cds_count}\t{i}\t{seq_record.name}\n")
                    seq_record.id = f"MG_{mg_count}-CDS_{cds_count}"
                    seq_record.name = ""
                    seq_record.description = ""
                    SeqIO.write(seq_record, fh, "fasta")


if config["params"]["dereplicate"]["cdhit"]["do_gene"]:
    rule dereplicate_gene_cdhit:
        input:
            ffn = os.path.join(config["output"]["predict"],
                               "scaftigs_gene_merged/{assembler}.prodigal.scaftigs.gene.merged.ffn"),
            metadata = os.path.join(config["output"]["predict"],
                                    "scaftigs_gene_merged/{assembler}.prodigal.scaftigs.gene.merged.ffn.metadata")
        output:
            os.path.join(config["output"]["dereplicate"],
                         "scaftigs_gene_nr/{assembler}.prodigal.scaftigs.gene.merged.nr.ffn.gz")
        conda:
            config["envs"]["cdhit"]
        log:
            os.path.join(config["output"]["dereplicate"],
                         "logs/cdhit/{assembler}.prodigal.scaftigs.gene.cdhit.log")
        benchmark:
            os.path.join(config["output"]["dereplicate"],
                         "benchmark/cdhit/{assembler}.cdhit.benchmark.txt")
        threads:
            config["params"]["dereplicate"]["cdhit"]["threads"]
        params:
            sequence_identity_threshold = config["params"]["dereplicate"]["cdhit"]["sequence_identity_threshold"],
            alignment_coverage_for_shorter_sequence = config["params"]["dereplicate"]["cdhit"]["alignment_coverage_for_shorter_sequence"],
            word_length = config["params"]["dereplicate"]["cdhit"]["word_length"],
            use_global_sequence_identity = config["params"]["dereplicate"]["cdhit"]["use_global_sequence_identity"],
            memory_limit = config["params"]["dereplicate"]["cdhit"]["memory_limit"],
            cluster_description_length = config["params"]["dereplicate"]["cdhit"]["cluster_description_length"],
            default_algorithm = config["params"]["dereplicate"]["cdhit"]["default_algorithm"],
            both_alignment = config["params"]["dereplicate"]["cdhit"]["both_alignment"]
        shell:
            '''
            FFN={output}
            ffn=`echo ${{FFN%.gz}}`

            cd-hit-est -i {input.ffn} -o $ffn \
            -c {params.sequence_identity_threshold} \
            -n {params.word_length} \
            -G {params.use_global_sequence_identity} \
            -aS {params.alignment_coverage_for_shorter_sequence} \
            -M {params.memory_limit} \
            -d {params.cluster_description_length} \
            -g {params.default_algorithm} \
            -r {params.both_alignment} \
            -T {threads} >{log} 2>&1

            pigz -p {threads} $ffn
            '''


    rule dereplicate_gene_cdhit_all:
        input:
            expand(os.path.join(
                config["output"]["dereplicate"],
                "scaftigs_gene_nr/{assembler}.prodigal.scaftigs.gene.merged.nr.ffn.gz"),
                   assembler=ASSEMBLERS)

else:
    rule dereplicate_gene_cdhit_all:
        input:


rule dereplicate_gene_all:
    input:
        rules.dereplicate_gene_cdhit_all.input
