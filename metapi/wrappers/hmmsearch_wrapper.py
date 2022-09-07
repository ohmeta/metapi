#!/usr/bin/env python
import pyhmmer
import sys

hmm_threads = sys.argv[0]
hmm_evalue = sys.argv[1]
hmm_tbl = sys.argv[2]
hmm_db = sys.argv[3]
hmm_seq = sys.argv[4]

# reference
# https://github.com/althonos/pyhmmer/issues/22

alphabet = pyhmmer.easel.Alphabet.amino()

with pyhmmer.easel.SequenceFile(hmm_seq, digital=True, alphabet=alphabet) as seq_file:
    sequences = list(seq_file)

with open(hmm_tbl, "wb") as dst:
    with pyhmmer.plan7.HMMFile(hmm_db) as hmm_file:
        for i, hits in enumerate(pyhmmer.hmmsearch(hmm_file, sequences, cpus=hmm_threads, E=hmm_evalue)):
            hits.write(dst, format="targets", header=i==0)

# example
# python hmmsearch_wrapper.py 8 0.01 output.tbl virus.hmm test.faa
