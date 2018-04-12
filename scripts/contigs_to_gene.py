#!/usr/bin/env python
#from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter
import shutil
import os
import errno
import argparse
import subprocess
import sys
import gzip

__author__ = 'Jie Zhu'
__version__ = '0.1.0'
__date__ = 'March 21, 2018'

def cut_fasta_by_len(fa_file, len_cutoff, outdir, prefix, suffix):
    # https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
    # Defeats race condition when another thread created the path
    #if not os.path.exists(outdir):
    #    os.mkdir(outdir)
    try:
        os.makedirs(outdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    cut_fa_file = os.path.join(outdir,
                               prefix + ".ge" + str(len_cutoff) + suffix)
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
def gene_prediction(cut_scatig_file, outdir, prefix, gz_or):
    '''gene prediction using prodigal
    
    https://github.com/hyattpd/Prodigal/wiki
    https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#metagenomes

    PRODIGAL v2.6.1 [July, 2013]
    Univ of Tenn / Oak Ridge National Lab
    Doug Hyatt, Loren Hauser, et al.
    -------------------------------------

    Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]

            -a:  Write protein translations to the selected file.
            -c:  Closed ends.  Do not allow genes to run off edges.
            -d:  Write nucleotide sequences of genes to the selected file.
            -f:  Select output format (gbk, gff, or sco).  Default is gbk.
            -g:  Specify a translation table to use (default 11).
            -h:  Print help menu and exit.
            -i:  Specify input file (default reads from stdin).
            -m:  Treat runs of n's as masked sequence and do not build genes across
                 them.
            -n:  Bypass the Shine-Dalgarno trainer and force the program to scan
                 for motifs.
            -o:  Specify output file (default writes to stdout).
            -p:  Select procedure (single or meta).  Default is single.
            -q:  Run quietly (suppress normal stderr output).
            -s:  Write all potential genes (with scores) to the selected file.
            -t:  Write a training file (if none exists); otherwise, read and use
                 the specified training file.
            -v:  Print version number and exit.
    '''

    pep_file = os.path.join(outdir, prefix + ".pep.faa")
    cds_file = os.path.join(outdir, prefix + ".cds.ffn")
    gff_file = os.path.join(outdir, prefix + ".cds.gff")
    start_file = os.path.join(outdir, prefix + ".score.gff")
    shell_file = os.path.join(outdir, prefix + ".prodigal.sh")

    prodigal = shutil.which("prodigal")
    shell = prodigal + \
            " -i " + cut_scatig_file + \
            " -a " + pep_file + \
            " -d " + cds_file + \
            " -o " + gff_file + \
            " -s " + start_file + \
            " -f gff -p meta -q"
    with open(shell_file, 'w') as out_f:
        out_f.write(shell + "\n")

    try:
        code = subprocess.check_output(shell, shell=True)
        print("gene prediction done!")
    except subprocess.CalledProcessError as e:
        code = e.returncode
        return code

    gzip_bin = shutil.which("gzip")
    if gz_or:
        gz_shell = gzip_bin + " " + cds_file + " && " + \
                   gzip_bin + " " + pep_file + " && " + \
                   gzip_bin + " " + gff_file + " && " + \
                   gzip_bin + " " + start_file
        try:
            subprocess.check_output(gz_shell, shell=True)
            print("gzip files done!")
        except subprocess.CalledProcessError as e:
            code = e.returncode
            return code

    return code

def main():
    parser = argparse.ArgumentParser(
        description="filter scaffolds or contigs by length cutoff and do gene prediction")
    parser.add_argument('-fa', type=str, help='scaffolds or contigs file path')
    parser.add_argument(
        '-sclen',
        type=int,
        help='scaffold or contigs length cutoff, default: 500',
        default=500)
    parser.add_argument(
        '-cdslen',
        type=int,
        help='cds length cutoff, default: 100',
        default=100)
    parser.add_argument(
        '-outdir', type=str, help='output dir store gene prediction results')
    parser.add_argument('-prefix', type=str, help='prefix for file name')
    parser.add_argument(
        '-rms',
        action='store_true',
        help=
        'delete cutoff scaffolds or contigs file after gene prediction, default: False, just -rm, it will true',
        default=False)
    parser.add_argument(
        '-gz',
        action='store_true',
        help='compress all output file, default: False, just -gz, it will true',
        default=False)
    args = parser.parse_args()

    # cut scaffolds or contigs
    cut_scatig_file = cut_fasta_by_len(args.fa, args.sclen, args.outdir, args.prefix + ".scatig", ".fa")
    
    # do gene prediction
    code = gene_prediction(cut_scatig_file, args.outdir, args.prefix, args.gz)

    # remove cutted scaffolds or contigs after gene prediction
    if code and args.rm:
        os.remove(cut_scatig_file)
        print("remove " + cut_scatig_file + " done!")

    # cut cds
    cds_file = os.path.join(args.outdir, args.prefix + ".cds.ffn")
    cut_cds_file = cut_fasta_by_len(cds_file, args.cdslen, args.outdir, args.prefix + ".cds", ".ffn")
    if cut_cds_file and args.gz:
        gz_shell = shutil.which("gzip") + " " + cut_cds_file 
        try:
            subprocess.check_output(gz_shell, shell=True)
        except subprocess.CalledProcessError as e:
            code = e.returncode
    
    # I think maybe we need another script to do this
    # we could cut gene and generate file by using gene file and gff file
    # finally, we get
    # cutted gene file
    # cutted gff file(maybe)

if __name__ == '__main__':
    main()