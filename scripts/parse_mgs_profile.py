#!/usr/bin/env python

import re
import sys
from pprint import pprint


def parse(mgs_profile):
    count = 0
    with open(mgs_profile, 'r') as ih:
        for line in ih:
            line_list = re.split(r"\s+|,", line)
            cag_id = line_list[0]
            seq_count = line_list[1]
            seq_id_list = line_list[2:]
            count += 1
            a = 0
            if count == 1:
                for i in seq_id_list:
                    print(i)
                    a += 1
                print(a)
                print(seq_count)
                break


def main():
    parse(sys.argv[1])


if __name__ == '__main__':
    main()
