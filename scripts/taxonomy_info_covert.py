#!/usr/bin/env python
import argparse
import csv
import os

def parse_lca_classify(taxonomy_csv, output):
    # taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']
    headers = ["ID", "status", "lineage"]
    rows = []
    with open(taxonomy_csv, 'r') as csv_h:
        f_csv = csv.DictReader(csv_h)
        # print(type(f_csv))
        for row in f_csv:
            row_dict = {}
            row_dict["ID"] = os.path.basename(row["ID"])
            row_dict["status"] = row["status"]
            row_dict["lineage"] = row["superkingdom"] + ";" + \
                                  row["phylum"] + ";" + \
                                  row["order"] + ";" + \
                                  row["class"] + ";" + \
                                  row["family"] + ";" + \
                                  row["genus"] + ";" + \
                                  row["species"]
            rows.append(row_dict)

    with open(output, 'w') as csv_out:
        csv_f = csv.DictWriter(csv_out, headers)
        csv_f.writeheader()
        csv_f.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description="convert sourmash lca classify resuts to metacoder input")
    parser.add_argument('-csv', type=str, help="sourmash lca classify results csv file")
    parser.add_argument('-out', type=str, help='coverted csv file')
    args = parser.parse_args()
    parse_lca_classify(args.csv, args.out)

if __name__ == '__main__':
    main()
