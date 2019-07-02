#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import time
import pickle
import subprocess
from Bio import bgzf
from taxadb.taxid import TaxID
from taxadb.names import SciName


def parse_kraken2_output(kraken2_output, prefix):
    start_time = time.time()
    taxid_counter = {}
    demultiplexer = {}

    with open(kraken2_output, 'r') as kh:
        for line in kh:
            cols = line.split('\t')
            read_id = cols[1]
            tax_name = cols[2].split("(")[0].strip()
            tax_id = int(cols[2].split("(")[-1].split(")")[0].split()[-1])

            demultiplexer[read_id] = tax_id
            if tax_id in taxid_counter:
                taxid_counter[tax_id][1] += 1
            else:
                taxid_counter[tax_id] = [tax_name, 1]

    print("step_1: parse kraken2 output has spent %d s" % (time.time() - start_time))

    with open(prefix.rstrip("/") + ".taxa.pickle", 'wb') as ph:
        pickle.dump(taxid_counter, ph)

    return taxid_counter, demultiplexer


def demultiplexer_kraken2(kraken2_output, r1, r2, change_seq_id, prefix, taxid_counter, demultiplexer):
    start_time = time.time()
    gzip_h = {}
    for i in taxid_counter:
        gzip_h[i] = {}
        gzip_h[i]["r1"] = bgzf.BgzfWriter(prefix.rstrip("/") + ".%d.1.fq.gz" % i, 'wb')
        gzip_h[i]["r2"] = bgzf.BgzfWriter(prefix.rstrip("/") + ".%d.2.fq.gz" % i, 'wb')

    if r1.endswith(".gz"):
        r1_h = gzip.open(r1, 'rt')
        r2_h = gzip.open(r2, 'rt')
    else:
        r1_h = open(r1, 'rt')
        r2_h = open(r2, 'rt')

    if change_seq_id:
        sample_tag = os.path.basename(prefix)

    for read_1, read_2 in zip(r1_h, r2_h):
        read_id = read_1[1:].split("/")[0]
        if change_seq_id:
            gzip_h[demultiplexer[read_id]]["r1"].write(">%s|%s%s%s%s" %
                (
                    sample_tag, read_1[1:], next(r1_h), next(r1_h), next(r1_h)
                ))
            gzip_h[demultiplexer[read_id]]["r2"].write(">%s|%s%s%s%s" %
                (
                    sample_tag, read_2[1:], next(r2_h), next(r2_h), next(r2_h)
                ))
        else:
            gzip_h[demultiplexer[read_id]]["r1"].write("%s%s%s%s" %
                (
                    read_1, next(r1_h), next(r1_h), next(r1_h)
                ))
            gzip_h[demultiplexer[read_id]]["r2"].write("%s%s%s%s" %
                (
                    read_2, next(r2_h), next(r2_h), next(r2_h)
                ))
    print("step_2: demultiplex taxid-reads has spent %d s" % (time.time() - start_time))

    for i in taxid_counter:
        gzip_h[i]["r1"].close()
        gzip_h[i]["r2"].close()
    r1_h.close()
    r2_h.close()


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
        '--prefix',
        default="./results/result",
        help='output prefix, default: ./'
    )
    parser.add_argument(
        '--rank',
        choices=["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
        default="genus",
        help='mini rank for merge'
    )
    parser.add_argument(
        '--change_seq_id',
        default=False,
        action='store_true',
        help='change seq id: add sample id to the front of seq id, default: False'
    )
    parser.add_argument(
        '--taxadb',
        type=str,
        help='taxonomy database'
    )
    parser.add_argument(
        '--merge_script',
        help='write merge command to merge script'
    )
    args = parser.parse_args(args_)

    LINEAGES = ["no_rank", "subspecies", "species", "genus", "family",
                "order", "class", "phylum", "superkingdom"]

    RANK = args.rank
    if not args.rank in LINEAGES[1:]:
        print("wrong rank %s" % args.rank)
        sys.exit(1)

    SUB_LINRAGES = LINEAGES[LINEAGES.index(RANK):]

    TAXID_DB = TaxID(dbtype='sqlite', dbname=args.taxadb)
    NAMES_DB = SciName(dbtype='sqlite', dbname=args.taxadb)

    os.makedirs(os.path.dirname(os.path.abspath(args.prefix)), exist_ok=True)

    def get_parent_taxid(tax_id, tax_name, level):
        if tax_id == 0:
            return "no_rank", 0, "unclassified"
    
        lineage_dict = TAXID_DB.lineage_id(tax_id, ranks=True)

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

    # step 1
    taxid_counter, demultiplexer = parse_kraken2_output(args.kraken2_output, args.prefix)

    # step 2
    start_time = time.time()
    ranks_counter = {}
    merger = {}
    total_reads_pair = 0
    for taxid in taxid_counter:
        rank, tax_id, taxa_name = get_parent_taxid(taxid, taxid_counter[taxid][0], RANK)
        if rank in ranks_counter:
            if tax_id in ranks_counter[rank]:
                ranks_counter[rank][tax_id][1] += taxid_counter[taxid][1]
            else:
                ranks_counter[rank][tax_id] = [taxa_name, taxid_counter[taxid][1]]
        else:
            ranks_counter[rank] = {}
            ranks_counter[rank][tax_id] = [taxa_name, taxid_counter[taxid][1]]

        if tax_id in merger:
            merger[tax_id]["r1"].append(os.path.abspath("%s.%d.1.fq.gz" % (args.prefix.rstrip("/"), taxid)))
            merger[tax_id]["r2"].append(os.path.abspath("%s.%d.2.fq.gz" % (args.prefix.rstrip("/"), taxid)))
        else:
            merger[tax_id] = {}
            merger[tax_id]["r1"] = [os.path.abspath("%s.%d.1.fq.gz" % (args.prefix.rstrip("/"), taxid))]
            merger[tax_id]["r2"] = [os.path.abspath("%s.%d.2.fq.gz" % (args.prefix.rstrip("/"), taxid))]

    print("step_2: summary kraken2 output has spent %d s" % (time.time() - start_time))
    print("rank\ttax_id\ttaxa_name\treads_count")
    for rank in LINEAGES:
        if rank in ranks_counter:
            print(rank)
            for tax_id in ranks_counter[rank]:
                print("\t%d\t%s\t%d" % (tax_id, ranks_counter[rank][tax_id][0], ranks_counter[rank][tax_id][1]))
                total_reads_pair += ranks_counter[rank][tax_id][1]
    print("total %d taxnonmy level" % len(ranks_counter))
    print("total %d reads pair" % total_reads_pair)
    
    # step 3
    #demultiplexer_kraken2(args.kraken2_output, args.r1, args.r2, args.change_seq_id, args.prefix, taxid_counter, demultiplexer)

    # step 4
    start_time = time.time()
    with open(args.merge_script, 'w') as oh:
        for tax_id in merger:
            r1_str = " ".join(merger[tax_id]["r1"])
            r2_str = " ".join(merger[tax_id]["r2"])
            
            if len(merger[tax_id]["r1"]) > 1:     
                cmd_1 = "cat %s > %s\n" % (r1_str, "%s.%d.merged.1.fq.gz" % (args.prefix.rstrip("/"), tax_id))
                cmd_2 = "cat %s > %s\n" % (r2_str, "%s.%d.merged.2.fq.gz" % (args.prefix.rstrip("/"), tax_id))
            else:
                cmd_1 = "mv %s %s\n" % (r1_str, "%s.%d.merged.1.fq.gz" % (args.prefix.rstrip("/"), tax_id))
                cmd_2 = "mv %s %s\n" % (r2_str, "%s.%d.merged.2.fq.gz" % (args.prefix.rstrip("/"), tax_id))
            
            rm_cmd_1 = "rm -rf %s\n" % r1_str
            rm_cmd_2 = "rm -rf %s\n" % r2_str
            
            oh.write(cmd_1)
            oh.write(rm_cmd_1)
            oh.write(cmd_2)
            oh.write(rm_cmd_2)
            #subprocess.call(cmd_1, shell=True)
            #subprocess.call(rm_cmd_1, shell=True)
            #subprocess.call(cmd_2, shell=True)
            #subprocess.call(rm_cmd_2, shell=True)
            
    print("step_4: merge taxid-reads has spent %d s" % (time.time() - start_time))


if __name__ == '__main__':
    main(sys.argv[1:])
