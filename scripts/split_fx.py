#!/usr/bin/env python
# please see http://biopython.org/wiki/Split_large_file
import argparse
import os
import errno

from Bio import SeqIO


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator. Each list will have
    batch_size entries, although the final list may be shorter. 
    """

    entry = True
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                # entry = iterator.next()
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                break
            batch.append(entry)
        if batch:
            yield batch


# TODO
# def split_fastq()
# def split_alignment() 


def split_fasta(fa_file, batch_size, outdir, onedir):
    record_iter = SeqIO.parse(open(fa_file, 'r'), "fasta")
    i = 0
    for i, batch in enumerate(batch_iterator(record_iter, batch_size), start = 1):
        if onedir:
            splitfa = os.path.join(outdir, "split_%i.fa" % (i))
        else:
            splitdir = os.path.join(outdir, "split_" + str(i))
            try:
                os.makedirs(splitdir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
            splitfa = os.path.join(splitdir, "split_%i.fa" % (i))
            
        with open(splitfa, 'w') as out_h:
            count = SeqIO.write(batch, out_h, "fasta")
        print("wrote %i records to %s" % (count, splitfa))
    return i


def main():
    """split large fasta/fastq file by seq size"""
    parser = argparse.ArgumentParser(description='split large fasta/fastq file by seq size')
    parser.add_argument('-f', type=str, help='input file, a large fasta or fastq file')
    parser.add_argument('-n', type=int, help='each splited file base size', default=1000)
    parser.add_argument('-outdir', type=str, help='a directory store splited file')

    args = parser.parse_args()
    split_fasta(args.f, args.n, args.outdir, False) 

if __name__ == '__main__':
    main()