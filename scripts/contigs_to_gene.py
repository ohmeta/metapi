#!/usr/bin/env python
#from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter
import shutil
import os
import argparse
import subprocess
import sys
import gzip


def filter_contigs_by_len(fa_file, len_cutoff, outdir, prefix):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    cut_fa_file = os.path.join(outdir,
                               prefix + "_gt" + str(len_cutoff) + ".fa")
    if os.path.exists(cut_fa_file) and (os.path.getsize(cut_fa_file) > 0):
        return cut_fa_file
    
    if fa_file.endswith(".gz"):
        in_h = gzip.open(fa_file, 'rt')
    else:
        in_h = open(fa_file, 'r')
    with open(cut_fa_file, 'w') as out_h:
        #for rec in SeqIO.parse(in_h, 'fasta'):
        #    if len(rec.seq) >= len_cutoff:
        #        SeqIO.write(rec, out_h, 'fasta')
        # yes, the SeqIO.parse() API is more simple to use, easy to understand
        # but, try different method, you will find something
        writer = FastaWriter(out_h)
        writer.write_header()
        for rec in FastaIterator(in_h):
            if len(rec) >= len_cutoff:
                writer.write_record(rec)
        writer.write_footer()
    in_h.close()
    return cut_fa_file


# just transfer seq(store in the RAM) to prodigal, the prodigal program also write seq to disk as tmp file
# so we store seq to disk first, then pass the seq file to prodigal by using -i parameter
def gene_prediction(cut_fa_file, len_cutoff, outdir, prefix, gz_or):
    '''gene prediction'''
    protein_file = os.path.join(outdir, prefix + ".protein.faa")
    cds_file = os.path.join(outdir, prefix + ".cds.ffn")
    cds_gff_file = os.path.join(outdir, prefix + ".cds.gff")
    start_file = os.path.join(outdir, prefix + ".score.gff")
    shell_file = os.path.join(outdir, prefix + ".prodigal.sh")

    prodigal = shutil.which("prodigal")
    shell = prodigal + \
            " -i " + cut_fa_file + \
            " -a " + protein_file + \
            " -d " + cds_file + \
            " -o " + cds_gff_file + \
            " -s " + start_file + \
            " -f gff -p anon -q"
    with open(shell_file, 'w') as out_f:
        out_f.write(shell + "\n")

    try:
        code = subprocess.check_output(shell, shell=True)
    except subprocess.CalledProcessError as e:
        code = e.returncode
        return code

    gzip_bin = shutil.which("gzip")
    if gz_or:
        gz_shell = gzip_bin + " " + protein_file + " && " + \
                   gzip_bin + " " + cds_file + " && " + \
                   gzip_bin + " " + cds_gff_file + " && " + \
                   gzip_bin + " " + start_file
        try:
            subprocess.check_output(gz_shell, shell=True)
        except subprocess.CalledProcessError as e:
            code = e.returncode
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="filter contigs and do gene prediction")
    parser.add_argument('-fa', type=str, help='contigs file path')
    parser.add_argument(
        '-len',
        type=int,
        help='contigs length cutoff, default: 500',
        default=500)
    parser.add_argument(
        '-outdir', type=str, help='output dir store gene prediction results')
    parser.add_argument('-prefix', type=str, help='prefix for file name')
    parser.add_argument(
        '-rm',
        action='store_true',
        help=
        'delete cutoff contigs after gene prediction, default: False, just -rm, it will true',
        default=False)
    parser.add_argument(
        '-gz',
        action='store_true',
        help='compress output file, default: False, just -gz, it will true',
        default=False)
    args = parser.parse_args()

    cut_fa_file = filter_contigs_by_len(args.fa, args.len, args.outdir,
                                        args.prefix)
    code = gene_prediction(cut_fa_file, args.len, args.outdir,
                           args.prefix, args.gz)

    if code == 0 and args.rm: 
        os.remove(cut_fa_file)


if __name__ == '__main__':
    main()