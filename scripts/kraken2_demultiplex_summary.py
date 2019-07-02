#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import pickle
from pprint import pprint
from taxadb.taxid import TaxID
from taxadb.names import SciName


def main(args_):
    parser = argparse.ArgumentParser("a summary of kraken2 demultiplex pickle")
    parser.add_argument(
        '--rank',
        choices=["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
        default="genus",
        help='mini rank for merge'
    )
    parser.add_argument(
        '--taxadb',
        type=str,
        help='taxonomy database'
    )
    parser.add_argument(
        '-p',
        '--pickle_list',
        help='kraken2 demultiplex pickle list'
    )
    parser.add_argument(
        '-o',
        '--summary_output'
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
        
    
    summary_dict = {"taxid": [],
                    "taxa_name": [],
                    "reads_count": [],
                    "rank": [],
                    "parent_taxid": [],
                    "parent_taxa_name": []}

    with open(args.pickle_list, 'r') as ih:
        for line in ih:
            with open(line.strip(), 'rb') as ph:
                kk2_ranks_counter = pickle.load(ph)
                # pprint(kk2_ranks_counter)

                for taxid in kk2_ranks_counter:
                    if taxid in summary_dict["taxid"]:
                        summary_dict["reads_count"][summary_dict["taxid"].index(taxid)] += 2 * kk2_ranks_counter[taxid][1]
                    else:
                        summary_dict["taxid"].append(taxid)
                        summary_dict["taxa_name"].append(kk2_ranks_counter[taxid][0])
                        summary_dict["reads_count"].append(2 * kk2_ranks_counter[taxid][1]) 
                        
                        rank_, taxid_, taxaname_ = get_parent_taxid(taxid, kk2_ranks_counter[taxid][0], RANK)
                        summary_dict["rank"].append(rank_)
                        summary_dict["parent_taxid"].append(taxid_)
                        summary_dict["parent_taxa_name"].append(taxaname_)

    summary_df = pd.DataFrame.from_dict(summary_dict)

    summary_df.to_csv(args.summary_output, index=False, sep='\t')


if __name__ == '__main__':
    main(sys.argv[1:])
