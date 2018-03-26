#!/usr/bin/env python
from Bio import SeqIO, bgzf
from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter
import gzip
import sys
import os
import argparse

#with open(sys.argv[2], 'w') as fa_out:
#    with open(sys.argv[1], 'r') as fa_in:
#        for rec in SeqIO.parse(fa_in, 'fasta'):
#            (description, sample_name) = rec.description.split("\t")
#            rec.description = sample_name + "_" + description
#            rec.id = rec.description.split(' ')[0]
#            SeqIO.write(rec, fa_out, 'fasta')

def change_header_sample(title):
    # title(total header) -> (id, name, description)
    # R0170300050_tooth_RA.contigs.fa
    #from > k119_1 flag=1 multi=7.0000 len=3284 R0170300050_tooth_RA
    #to   > R0170300050_tooth_RA_k119_1 flag=1 multi=7.0000 len=3284
    (one_line, sample_name) = title.split("\t")
    id = sample_name + "_" + title.split(' ')[0]
    desc = id + ' ' + ' '.join(one_line.split(' ')[1:])
    return id, "", desc

def change_header_no_sample(title):
    # title(total header) -> (id, name, description)
    # R0170300050_tooth_RA.contigs.fa
    #from > k119_1 flag=1 multi=7.0000 len=3284
    #to   > R0170300050_tooth_RA_k119_1 flag=1 multi=7.0000 len=3284
    id = sample_tag + "_" + title.split(' ')[0]
    desc = id + " " + " ".join(title.split(' ')[1:])
    return id, "", desc


## rename header framework
## just change header_function
def reheader_fasta(fa_in, fa_out, header_function, in_gz, gz):
    if in_gz:
        in_h = gzip.open(fa_in, 'rt')
    else:
        in_h = open(fa_in, 'r')
    if gz:
        out_h = bgzf.BgzfWriter(fa_out, 'wb')
    else:
        out_h = open(fa_out, 'w')
    writer = FastaWriter(out_h)
    writer.write_header()
    for rec in FastaIterator(in_h, title2ids = header_function):
        writer.write_record(rec)
    writer.write_footer()
    out_h.close()
    in_h.close()

def main():
    '''
    Why write this script ?
    Becaust megahit always generate knum_num format contigs id
    '''
    parser = argparse.ArgumentParser(description='change fasta file header')
    parser.add_argument('-fa', type=str, help='fasta file path')
    parser.add_argument('-out', type=str, help='output')
    parser.add_argument('-rm', action='store_true', help='delete original fasta file', default=False)
    parser.add_argument('-gz', action='store_true', help='compress output fasta file', default=False)
    parser.add_argument('-mv', action='store_true', help="rename change id fasta file to original file", default=False)
    args = parser.parse_args()
    
    #assert not args.fa == args.out, "input file name can't equal to output file name"
    if (args.out == args.fa) or (not args.out):
        args.out = args.fa + ".changeid"
    if args.gz:
        if not args.out.endswith(".gz"):
            args.out = args.out + ".gz"
    
    in_gz = args.fa.endswith(".gz")
    #if args.fa.endswith(".gz"):
    #    args.gz = True
    
    global sample_tag
    sample_tag = os.path.basename(args.fa).split(".")[0]
    
    abs_in = os.path.abspath(args.fa)
    abs_out = os.path.abspath(args.out)
    reheader_fasta(abs_in, abs_out, change_header_no_sample, in_gz, args.gz)
    
    if args.rm:
        os.remove(abs_in)
    if args.mv:
        if (not in_gz) and args.gz:
            abs_in = abs_in + ".gz"
        if in_gz and (not args.gz):
            abs_in = abs_in.rstrip(".gz")
        os.rename(abs_out, abs_in)

if __name__ == '__main__':
    main()