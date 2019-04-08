#!/usr/bin/env python
import pandas as pd
import sys
import os


def main():
    samples = pd.read_csv(sys.argv[1], sep='\t').set_index("id", drop=False)
    for i in samples.index:
        fq1, fq2 = samples.loc[i, ["fq1", "fq2"]]
        if (not os.path.exists(fq1)) or (not os.path.exists(fq2)):
            print("error: %s\t%s\t%s" % (i, fq1, fq2))


if __name__ == '__main__':
    main()
