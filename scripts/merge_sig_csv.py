#!/usr/bin/env python

import csv
import argparse
import pandas as pd

def merge_csv(csvlist, output):
    dfs = pd.DataFrame([])
    with open(csvlist, 'r') as csv_l:
        for csv_f in csv_l:
            dfs.append(pd.read_csv(csv_f.strip()))
    dfs.to_csv(output)

def main():
    parser = argparse.ArgumentParser(description="merge sourmash sigs to a csv file")
    parser.add_argument('-csvlist', type=str, help='a file contain sig file path list')
    parser.add_argument('-output', type=str, help='output csv file')
    args = parser.parse_args()

if __name__ == '__main__':
    main()
