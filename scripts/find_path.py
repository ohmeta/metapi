#!/usr/bin/env python3
import os

def find_path(dir, suffix):
    path = {}
    for f in os.listdir(dir):
        if f.endswith(suffix):
            key = f.rstrip(suffix)
            path[key] = os.path.join(dir, f)
    return path

def find_path_tag(dir, tag):
    if tag == "raw":
        r1 = {}
        r2 = {}
        for f in os.listdir(dir):
            if f.endswith("1.fq.gz"):
                key = f.rstrip(".|-|_" + "1.fq.gz")
                r1[key] = os.path.join(dir, f)
            if f.endswith("2.fq.gz"):
                key = f.rstrip(".|-|_" + "2.fq.gz")
                r2[key] = os.path.join(dir, f)
        return (r1, r2)
    elif tag == "clean" or tag == "rmhost":
        r1 = {}
        r2 = {}
        rs = {}
        rt = {}
        for f in os.listdir(dir):
            if f.endswith(tag + ".1.fq.gz"):
                key = f.rstrip(".|-|_" + tag + ".1.fq.gz")
                r1[key] = os.path.join(dir, f)
            if f.endswith(tag + ".2.fq.gz"):
                key = f.rstrip(".|-|_" + tag + ".2.fq.gz")
                r1[key] = os.path.join(dir, f)
            if f.endswith(tag + ".single.fq.gz"):
                key = f.rstrip(".|-|_" + tag + ".single.fq.gz")
                r1[key] = os.path.join(dir, f)
            if f.endswith(tag + ".stat_out"):
                key = f.rstrip(".|-|_" + tag + ".stat_out")
                r1[key] = os.path.join(dir, f)
        return (r1, r2, rs, rt)
