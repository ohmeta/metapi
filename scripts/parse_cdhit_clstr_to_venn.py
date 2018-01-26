#!/usr/bin/env python
import argparse
import re
import sys


def parse_cdhit_clstr(clstr_f):
    '''parse cdhit cluster file'''
    geneset = {}
    cluster_num = 0
    id_overlap = []
    rat_n = 0
    mouse_n = 0
    human_n = 0

    rat_tag = False
    mouse_tag = False
    human_tag = False

    geneset["rat"] = 0   # total rat
    geneset["mouse"] = 0   # total mouse
    geneset["human"] = 0   # total human
    
    geneset["rat_mouse"] = 0   # only rat and mouse
    geneset["rat_human"] = 0   # only rat and human
    geneset["mouse_human"] = 0   # only mouse and human
    geneset["rat_mouse_human"] = 0   # only rat and mosue and human
    
    geneset["rat_mouse_rat"] = 0
    geneset["rat_mouse_mouse"] = 0
    geneset["rat_human_rat"] = 0
    geneset["rat_human_human"] = 0
    geneset["mouse_human_mouse"] = 0
    geneset["mouse_human_human"] = 0
    
    geneset["rat_mouse_human_rat"] = 0
    geneset["rat_mouse_human_mouse"] = 0
    geneset["rat_mouse_human_human"] = 0

    pattern = re.compile(r">(.+)\.\.\.")

    with open(clstr_f, 'r') as handle:
        for i in handle:
            if i.startswith(">"):
                cluster_num += 1
                if id_overlap:
                    for gene_id in id_overlap:
                        if "rat" in gene_id:
                            geneset["rat"] += 1
                            rat_n += 1
                            rat_tag = True
                        elif "mouse" in gene_id:
                            geneset["mouse"] += 1
                            mouse_n += 1
                            mouse_tag = True
                        elif "human" in gene_id:
                            geneset["human"] += 1
                            human_n += 1
                            human_tag = True
                        else:
                            sys.exit("wrong id")
                    
                    if (rat_tag and mouse_tag) and (not human_tag):
                        geneset["rat_mouse"] += 1
                        geneset["rat_mouse_rat"] += rat_n
                        geneset["rat_mouse_mouse"] += mouse_n
                    
                    if (rat_tag and human_tag) and (not mouse_tag):
                        geneset["rat_human"] += 1
                        geneset["rat_human_rat"] += rat_n
                        geneset["rat_human_human"] += human_n
                    
                    if (mouse_tag and human_tag) and (not rat_tag):
                        geneset["mouse_human"] += 1
                        geneset["mouse_human_mouse"] += mouse_n
                        geneset["mouse_human_human"] += human_n
                    
                    if rat_tag and mouse_tag and human_tag:
                        geneset["rat_mouse_human"] += 1
                        geneset["rat_mouse_human_rat"] += rat_n
                        geneset["rat_mouse_human_mouse"] += mouse_n
                        geneset["rat_mouse_human_human"] += human_n
               
                id_overlap = []
                rat_n = 0
                mouse_n = 0
                human_n = 0
                rat_tag = False
                mouse_tag = False
                human_tag = False
            else:
                matches = re.search(pattern, i)
                if matches:
                    id_overlap.append(matches.group(1))
                else:
                    sys.exit("can't extract id")

        if id_overlap:
            for gene_id in id_overlap:
                if "rat" in gene_id:
                    geneset["rat"] += 1
                    rat_n += 1
                    rat_tag = True
                elif "mouse" in gene_id:
                    geneset["mouse"] += 1
                    mouse_n += 1
                    mouse_tag = True
                elif "human" in gene_id:
                    geneset["human"] += 1
                    human_n += 1
                    human_tag = True
                else:
                    sys.exit("wrong id")
                    
            if (rat_tag and mouse_tag) and (not human_tag):
                geneset["rat_mouse"] += 1
                geneset["rat_mouse_rat"] += rat_n
                geneset["rat_mouse_mouse"] += mouse_n
                    
            if (rat_tag and human_tag) and (not mouse_tag):
                geneset["rat_human"] += 1
                geneset["rat_human_rat"] += rat_n
                geneset["rat_human_human"] += human_n
                    
            if (mouse_tag and human_tag) and (not rat_tag):
                geneset["mouse_human"] += 1
                geneset["mouse_human_mouse"] += mouse_n
                geneset["mouse_human_human"] += human_n
                    
            if rat_tag and mouse_tag and human_tag:
                geneset["rat_mouse_human"] += 1
                geneset["rat_mouse_human_rat"] += rat_n
                geneset["rat_mouse_human_mouse"] += mouse_n
                geneset["rat_mouse_human_human"] += human_n

    
    print("total cluster: %d\n" % cluster_num)
    
    print("total rat gene number: %d" % geneset["rat"])
    print("total mouse gene number: %d" % geneset["mouse"])
    print("total human gene number: %d\n" % geneset["human"])
    
    print("total rat and mouse cluster: %d" % geneset["rat_mouse"])
    print("\trat gene number: %d" % geneset["rat_mouse_rat"])
    print("\tmouse gene number: %d\n" % geneset["rat_mouse_mouse"])

    print("total rat and human cluster: %d" % geneset["rat_human"])
    print("\trat gene number: %d" % geneset["rat_human_rat"])
    print("\thuman gene number: %d\n" % geneset["rat_human_human"])

    print("total mouse and human cluster: %d" % geneset["mouse_human"])
    print("\tmouse gene number: %d" % geneset["mouse_human_mouse"])
    print("\thuman gene number: %d\n" % geneset["mouse_human_human"])

    print("total rat and mouse and human cluster: %d" % geneset["rat_mouse_human"])
    print("\trat gene number: %d" % geneset["rat_mouse_human_rat"])
    print("\tmouse gene number: %d" % geneset["rat_mouse_human_mouse"])
    print("\thuman gene number: %d\n" % geneset["rat_mouse_human_human"])

    return geneset


def main():
    '''main function'''
    parser = argparse.ArgumentParser(
        description="parse cdhit clstr and plot venn")
    parser.add_argument('-f', '--clstr', type=str,
                        help='cdhit cluster file')
    #parser.add_argument('-p', '--prefix', type=str,
    #                    help='output prefix')
    args = parser.parse_args()

    parse_cdhit_clstr(args.clstr)


if __name__ == '__main__':
    main()
