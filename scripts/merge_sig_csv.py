#!/usr/bin/env python

import csv
import argparse
import pandas as pd

def merge_csv(csvlist, output):
    frame = pd.DataFrame()
    list = []
    with open(csvlist, 'r') as csv_l:
        for csv_f in csv_l:
            df = pd.read_csv(csv_f.strip(), index_col=None, header=0)
            list.append(df)
    frame = pd.concat(list)
    frame.to_csv(output)

def main():
    parser = argparse.ArgumentParser(description="merge sourmash sigs to a csv file")
    parser.add_argument('-csvlist', type=str, help='a file contain sig file path list')
    parser.add_argument('-output', type=str, help='output csv file')
    args = parser.parse_args()
    merge_csv(args.csvlist, args.output)

if __name__ == '__main__':
    main()
