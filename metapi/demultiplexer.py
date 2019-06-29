#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import time
from datetime import datetime
import pandas as pd
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from taxadb.taxid import TaxID
from taxadb.names import SciName


def main(args_):
    parser = argparse.ArgumentParser("demultiplex PE fastq by paired reads classification of kraken2 output")
    parser.add_argument(
        '--r1',
        help='r1 fastq input'
    )
    parser.add_argument(
        '--r2',
        help='r2 fastq input'
    )
    parser.add_argument(
        '--kraken2_output',
        help='kraken2 output file'
    )
    parser.add_argument(
        '--rank',
        default="genus",
        choices=["superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"],
        help='highest taxonomy level, default: genus'
    )
    parser.add_argument(
        '--taxadb',
        help='ncbi taxaonomy lineages sqlite database'
    )
    parser.add_argument(
        '--prefix',
        default="./",
        help='output prefix, default: ./'
    )
    parser.add_argument(
        '--change_seq_id',
        default=False,
        action='store_true',
        help='change seq id: add sample id to the front of seq id, default: False'
    )
    parser.add_argument(
        '--log',
        help='log output'
    )
    parser.add_argument(
        '--debug',
        default=False,
        action='store_true',
        help='debug it?'
    )
    args = parser.parse_args(args_)

    LINEAGES = ["no_rank", "subspecies", "species", "genus", "family",
                "order", "class", "phylum", "superkingdom"]

    RANK = "genus"
    if not args.rank in LINEAGES[1:]:
        print("wrong rank %s" % args.rank)
        sys.exit(1)
    elif args.rank is None:
        RANK = "genus"

    TAXID_DB = TaxID(dbtype='sqlite', dbname=args.taxadb)
    NAMES_DB = SciName(dbtype='sqlite', dbname=args.taxadb)

    SUB_LINRAGES = LINEAGES[LINEAGES.index(RANK):]

    def get_parent_taxid(tax_id, tax_name, level):
        if tax_id == 0:
            return "no_rank", 0, "unclassified"
    
        lineage_dict = TAXID_DB.lineage_id(tax_id, ranks=True)

        # In kraken2 output, we can meet 2071623
        
        # > TAXID_DB.lineage_id(2071623, ranks=True) is None
        # > True

        # > NAMES_DB.taxid("Laceyella sp. FBKL4.010") is None
        # > True

        # > NAMES_DB.taxid("Laceyella sp. FBKL4.010".split()[0])
        # 292635

        # TAXID_DB.lineage_id(292635, ranks=True)
        # > {'genus': 292635,
        #    'family': 186824,
        #    'order': 1385,
        #    'class': 91061,
        #    'phylum': 1239,
        #    'no rank': 131567,
        #    'superkingdom': 2}

        if lineage_dict is None:
            taxid = NAMES_DB.taxid(tax_name)
            if taxid is None:
                taxid = NAMES_DB.taxid(tax_name.split()[0])
            if not taxid is None:
                lineage_dict = TAXID_DB.lineage_id(taxid, ranks=True)
            else:
                return "no_rank", tax_id, tax_name
    
        for rank in SUB_LINRAGES:
            if rank in lineage_dict:
                return rank, lineage_dict[rank], TAXID_DB.lineage_name(lineage_dict[rank], ranks=True)[rank]
        return "no_rank", tax_id, "unclassified"


    def parse_taxonomy_lineage(lineage_csv, taxid="tax_id"):
        taxonomy_lineages = pd.read_csv(lineage_csv, sep='\t', low_memory=False)\
                              .set_index(taxid, drop=False)

    start_time = time.time()

    ranks_counter = {}
    demultiplexer = {}

    os.makedirs(os.path.dirname(os.path.abspath(args.prefix)), exist_ok=True)

    log_h = open(args.log, 'w')
    log_h.write("now: %s\n" % datetime.today())

    count = 0
    with open(args.kraken2_output, 'r') as kh:
        for line in kh:
            cols = line.split('\t')
            read_id = cols[1]
            tax_name = cols[2].split("(")[0].strip()
            tax_id = int(cols[2].split("(")[-1].split(")")[0].split()[-1])
 
            if args.debug:
                print("%s %s %d" % (read_id, tax_name, tax_id))
            
            rank, tax_id_, taxa_name_ = get_parent_taxid(tax_id, tax_name, RANK)
            
            if rank in ranks_counter:
                if tax_id_ in ranks_counter[rank]:
                    ranks_counter[rank][tax_id_][1] += 1
                else:
                    ranks_counter[rank][tax_id_] = [taxa_name_, 1]
            else:
                ranks_counter[rank] = {}
                ranks_counter[rank][tax_id_] = [taxa_name_, 1]

            if tax_id_ in demultiplexer:
                demultiplexer[tax_id_].append(read_id)
            else:
                demultiplexer[tax_id_] = [read_id]
            #if args.debug:
            #    count += 1
            #    if count >= 1000:
            #        #break
            #        sys.exit(1)

    log_h.write("step_1: parse kraken2 output has spent %d s\n" % (time.time() - start_time))
    
    total_reads_pair = 0
    log_h.write("\trank\ttaxid\ttaxaonomy\tcount\n")
    for rank in LINEAGES:
        if rank in ranks_counter:
            log_h.write("\t%s\n" % rank)
            for tax_id_ in ranks_counter[rank]:
                log_h.write("\t\t%d\t%s\t%d\n" % (tax_id_, 
                                                  ranks_counter[rank][tax_id_][0],
                                                  ranks_counter[rank][tax_id_][1]))
                total_reads_pair += ranks_counter[rank][tax_id_][1]

    log_h.write("total %d taxnonmy level\n" % len(ranks_counter))
    log_h.write("total %d reads pair\n" % total_reads_pair)
    s1_time = time.time()
    
    log_h.write("step_2: begin build index for paired reads\n")

    def set_key(key):
        return key.split("/")[0]

    sample_tag = os.path.basename(args.prefix)

    db1 = args.prefix + ".1.db"
    db2 = args.prefix + ".2.db"
    
    if os.path.exists(db1):
        os.remove(db1)
    if os.path.exists(db2):
        os.remove(db2)

    r1_db = SeqIO.index_db(db1, args.r1, "fastq", key_function=set_key)
    r2_db = SeqIO.index_db(db2, args.r2, "fastq", key_function=set_key)
    
    log_h.write("step_2: build index has spent %d s\n" % (time.time() - s1_time))
    log_h.write("total %d reads_1 and %d reads_2\n" % (len(r1_db), len(r2_db)))
    
    s2_time = time.time()

    log_h.write("step_3: begin demultiplex taxid-reads\n")
    
    for tax_id in demultiplexer:
        with bgzf.BgzfWriter(args.prefix + ".%d.1.fq.gz" % tax_id, 'wb') as r1_oh, \
             bgzf.BgzfWriter(args.prefix + ".%d.2.fq.gz" % tax_id, 'wb') as r2_oh:
            for read_id in demultiplexer[tax_id]:
                if args.change_seq_id:
                    '''
                    r1_record = SeqRecord(r1_db[read_id].seq, 
                                          id=sample_tag + "|" + r1_db[read_id].id,
                                          description=sample_tag + "|" + r1_db[read_id].description,
                                          letter_annotations=r1_db[read_id].letter_annotations)
                    r2_record = SeqRecord(r2_db[read_id].seq, 
                                          id=sample_tag + "|" + r2_db[read_id].id,
                                          description=sample_tag + "|" + r2_db[read_id].description,
                                          letter_annotations=r2_db[read_id].letter_annotations)
                    '''
                    r1_oh.write("@%s|%s\n%s\n+\n%s\n" % \
                        (sample_tag, r1_db[read_id].description,
                         str(r1_db[read_id].seq),
                         SeqIO.QualityIO._get_sanger_quality_str(r1_db[read_id])))
                    r2_oh.write("@%s|%s\n%s\n+\n%s\n" % \
                        (sample_tag,
                         r2_db[read_id].description,
                         str(r2_db[read_id].seq),
                         SeqIO.QualityIO._get_sanger_quality_str(r2_db[read_id])))
                else:
                    SeqIO.write(r1_db[read_id], r1_oh, 'fastq')
                    SeqIO.write(r2_db[read_id], r2_oh, 'fastq')

    log_h.write("demultiplex %d paired reads has spent %d s\n" % (len(r1_db), time.time() - s2_time))
    log_h.write("kraken2 total reads pair: %d\n" % total_reads_pair)
    log_h.write("fastq total reads pair: %d\n" % len(r1_db))
    log_h.write("%d equal %d ?: %s\n"  % (len(r1_db), total_reads_pair, "True" if len(r1_db) == total_reads_pair else "False"))
    log_h.write("total time %d s\n" % (time.time() - start_time))
    log_h.write("now: %s\n" % datetime.today())

    log_h.close()


if __name__ == '__main__':
    main(sys.argv[1:])
