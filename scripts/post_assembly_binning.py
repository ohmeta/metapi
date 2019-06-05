#!/usr/bin/env python
import argparse
import os

import pandas as pd


def codegen(samples_tsv, output_dir):
    samples = pd.read_csv(samples_tsv, sep='\t').set_index("bin_id", drop=False)
    os.makedirs(output_dir, exist_ok=True)
    index_dir = os.path.join(output_dir, "00.index")
    os.makedirs(index_dir, exist_ok=True)
    mapping_dir = os.path.join(output_dir, "01.mapping")
    os.makedirs(mapping_dir, exist_ok=True)
    asm_dir = os.path.join(output_dir, "02.assembly")
    os.makedirs(asm_dir, exist_ok=True)

    with open(os.path.join(output_dir, "step1.index.sh"), 'w') as oh1, \
        open(os.path.join(output_dir, "step2.mapping.sh"), 'w') as oh2, \
        open(os.path.join(output_dir, "step3.assembly_spades.sh"), 'w') as oh3, \
        open(os.path.join(output_dir, "step3.assembly_shovill_spades.sh"), 'w') as oh4, \
        open(os.path.join(output_dir, "step3.assembly_shovill_megahit.sh"), 'w') as oh5, \
        open(os.path.join(output_dir, "step3.assembly_shovill_velvet.sh"), 'w') as oh6, \
        open(os.path.join(output_dir, "step3.assembly_shovill_skesa.sh"), 'w') as oh7:
        for bin_id in samples.index:
            # index
            prefix = os.path.join(index_dir, bin_id)
            cmd = "bwa index %s -p %s\n" % (samples.loc[bin_id, "bins_fna_path"], prefix)
            oh1.write(cmd)

            # mapping and extract reads
            r1 = os.path.join(mapping_dir, "%s.r1.fq.gz" % bin_id)
            r2 = os.path.join(mapping_dir, "%s.r2.fq.gz" % bin_id)
            stat = os.path.join(mapping_dir, "%s-flagstat.txt" % bin_id)
            cmd = "bwa mem -t 8 %s %s %s | tee >(samtools flagstat -@8 - > %s) | samtools fastq -@8 -F 12 -n -1 %s -2 %s -\n" % (
                prefix, samples.loc[bin_id, "fq1"], samples.loc[bin_id, "fq2"], stat, r1, r2)
            oh2.write(cmd)

            # assembly
            ## spades
            bins_asm_dir = os.path.join(asm_dir, bin_id + ".spades_out")
            cmd = "spades.py -1 %s -2 %s -k 21,29,39,59,79,99 --threads 8 -o %s\n" % (r1, r2, bins_asm_dir)
            oh3.write(cmd)
            ## shovill
            ### spades or megahit or velvet or skesa
            for assembler, file_handle in zip(["spades", "megahit", "velvet", "skesa"], [oh4, oh5, oh6, oh7]):
                bins_asm_dir = os.path.join(asm_dir, bin_id + ".shovill_%s_out" % assembler)
                cmd = "shovill --cpus 8 --keepfiles --assembler %s --outdir %s --R1 %s --R2 %s\n" % (assembler, bins_asm_dir, r1, r2)
                file_handle.write(cmd)


def main():
    parser = argparse.ArgumentParser(description='reassembly reads')
    parser.add_argument('-s', '--samples', type=str, help='metagenomics bins and paired reads list')
    parser.add_argument('-o', '--outdir', type=str, help='output directory')
    args = parser.parse_args()
    codegen(args.samples, args.outdir)


if __name__ == "__main__":
    main()
