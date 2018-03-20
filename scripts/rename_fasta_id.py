#!/usr/bin/env python
from Bio import SeqIO
import sys

with open(sys.argv[2], 'w') as fa_out:
    with open(sys.argv[1], 'r') as fa_in:
        for rec in SeqIO.parse(fa_in, 'fasta'):
            #print("id:" + rec.id)
            #print("description: " + rec.description)
            (description, sample_name) = rec.description.split("\t")
            rec.description = sample_name + "_" + description
            # 改变description并不会影响到之前的id,所以要根据id的定义对id重新赋值
            rec.id = rec.description.split(' ')[0]
            #print("id:" + rec.id)
            #print("description: " + rec.description)
            #rec.id = rec.description.split(' ')[0]
            #print("id:" + rec.id)
            
            #print(rec)
            # name会随着description的改变而跟着变化，但是要等到再次解析这条seq的时候才起作用
            #print(rec.name)
            SeqIO.write(rec, fa_out, 'fasta')

#with open(sys.argv[2], 'r') as fa_in2:
#    for rec in SeqIO.parse(fa_in2, 'fasta'):
#        print(rec)

# 所以，整体上看，SeqIO.parse()函数注重解析，在改变序列的属性上，逻辑关系处理得不好

def change_id(title):
    # title(total header) -> (id, name, description)
    # > k119_1 flag=1 multi=7.0000 len=3284 R0170300050_tooth_RA
    # > R0170300050_tooth_RA_k119_1 flag=1 multi=7.0000 len=3284
    (one_line, sample_name) = title.split("\t")
    id = sample_name + "_" + title.split(' ')[0]
    desc = id + ' ' + ' '.join(one_line.split(' ')[1:])
    return id, "", desc


def rename_fasta(fasta_f, sample_name, out):
    from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter
    with open(fasta_f, 'r') as fa_handle, open(out, 'w') as out_h:
        fa_writer = FastaWriter(out_h)
        for rec in FastaIterator(fa_handle, title2ids = change_id):
            fa_writer.write_file(rec)