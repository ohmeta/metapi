#!/usr/bin/env python
# please see https://github.com/lh3/asub 

import argparse
import fileinput
import subprocess

def main():
    '''it is a very simple script to submit array job'''
    with fileinput.input() as f_input:
        for line in f_input:
            print(line, end='')



if __name__ == '__main__':
    main()