#!/usr/bin/env python
import sys


def contigs_from_sample(contigs_len, sc_out):
    info = {}
    #count = 0
    with open(contigs_len, 'r') as handle:
        for line in handle:
            key = '_'.join(line.split("_")[:3])
            len = int(line.split("\t")[-1])
            if key in info:
                info[key]["num"] += 1
                info[key]["len"] += len
            else:
                info[key] = {}
                info[key]["num"] = 0
                info[key]["len"] = 0
            #count += 1
            #if count == 10000:
            #    break
    with open(sc_out, 'w') as out:
        out.write("sample_name\ttotal_contigs_num\ttotal_contigs_len\n")
        for key in info:
            out.write(key + "\t" + str(info[key]["num"]) + "\t" +
                      str(info[key]["len"]) + "\n")


def main():
    contigs_from_sample(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
