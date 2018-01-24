#!/usr/bin/python
'''get .fq.fqStat.txt, .clean.stat_out, .rmhost.stat_out info'''
import argparse
import csv
import os
import re
#import glob
#import sys

import numpy as np

def parse_fqstat(handle):
    '''parse .fq.fqStat.txt'''
    line_dict = {}
    num_list = []
    for line in handle:
        line_list = line.strip().split('\t')
        if line.startswith('#'):
            key = line_list[0].lstrip('#')
            if key != "Pos":
                if key == "Name":
                    line_dict[key] = line_list[1].split("\\")[-1].split(".")[0]
                else:
                    line_dict[key] = line_list[1]
                    if key == "N_Count":
                        line_dict["N_Rate"] = line_list[2]
        else:
            #onerow_num = [int(num) for num in line_list[6:-1:1]
            #num_list.append(onerow_num)
            num_list.append(line_list[6:-1:1])
    base_qual_mat = np.array(num_list, dtype='int64')
    #print(base_qual_mat.shape)
    return (line_dict, base_qual_mat)

def merge_fqstat(reads_type, dir_name, reads_len=100, max_qual=42):
    '''merge fqstat'''
    path_list = []
    #os.chdir(dir_name)
    #for file in glob.glob("*.fqStat.txt"):
    #    path_list.append(file)
    for file in os.listdir(dir_name):
        if file.endswith("fqStat.txt"):
            path_list.append(os.path.join(dir_name, file))

    line_dict = {}
    fqstat_dict = {}
    fqstat_dict["headers"] = []
    fqstat_dict["bodyinfo"] = []
    qual_mat_suma = np.zeros((reads_len, max_qual), dtype='int64')
    qual_mat_sumb = np.zeros((reads_len, max_qual), dtype='int64')

    for fq_path in path_list:
        fq_path = fq_path.strip()
        with open(fq_path, 'r') as handle:
            fqstat_tuple = parse_fqstat(handle)
            line_dict = fqstat_tuple[0]
            #print(fqstat_tuple[1].shape)
            if reads_type == "SE" or (reads_type == "PE" and
                                      re.search(r"_1.fq.fqStat.txt$", fq_path)):
                qual_mat_suma += fqstat_tuple[1]
            else:
                qual_mat_sumb += fqstat_tuple[1]
            fqstat_dict["bodyinfo"].append(line_dict)
    fqstat_dict["headers"] = line_dict.keys()
    sorted(fqstat_dict["headers"])

    #script_dir = os.path.abspath(sys.path[0])
    #os.chdir(script_dir)

    if reads_type == "SE":
        return (fqstat_dict, qual_mat_suma)
    if reads_type == "PE":
        return (fqstat_dict, qual_mat_suma, qual_mat_sumb)

def write_csv(writefile, headers, bodyinfo):
    '''write dict to csv file'''
    with open(writefile, 'w') as handle:
        f_csv = csv.DictWriter(handle, headers)
        f_csv.writeheader()
        f_csv.writerows(bodyinfo)

def write_qual_mat(mat_out, qual_mat_sum):
    '''write matrix object to txt file'''
    np.savetxt(mat_out, qual_mat_sum)

def get_fqstat(reads_type, path_list, fqstat_out):
    '''get fqstat info'''
    dir_name = os.path.dirname(fqstat_out)
    file_name = os.path.splitext(os.path.basename(fqstat_out))[0]

    if reads_type == "SE":
        fqstat_tuple = merge_fqstat(reads_type, path_list)
        write_csv(fqstat_out, fqstat_tuple[0]["headers"], fqstat_tuple[0]["bodyinfo"])
        mat_name = os.path.join(dir_name, file_name + ".mat.txt")
        write_qual_mat(mat_name, fqstat_tuple[1])

    if reads_type == "PE":
        fqstat_tuple = merge_fqstat(reads_type, path_list, 100, 42)
        write_csv(fqstat_out, fqstat_tuple[0]["headers"], fqstat_tuple[0]["bodyinfo"])
        mat_namea = os.path.join(dir_name, file_name + "_mata.txt")
        mat_nameb = os.path.join(dir_name, file_name + "_matb.txt")
        write_qual_mat(mat_namea, fqstat_tuple[1])
        write_qual_mat(mat_nameb, fqstat_tuple[2])

def parse_statout(handle, statout_path, flag="clean", reads_type="PE"):
    '''parse .clean.stat_out or .rmhost.stat_out'''
    col_dict = {}
    if flag == "clean":
    	headers_list = [key for key in handle.readline().strip().split("\t")[0:-1]]
    else:
    	headers_list = [key for key in handle.readline().strip().split("\t")]
    bodyinfo_list = [value for value in handle.readline().strip().split("\t")]
    
    headers_list.append("Name")
    statout_name = os.path.basename(statout_path).split(".", maxsplit = 1)[0]
    bodyinfo_list.append(statout_name)
    
    if flag == "clean" and reads_type == "PE":
        for key in handle.readline().strip().split("\t"):
            headers_list.append(key + "_1")
        for val in handle.readline().strip().split("\t"):
            bodyinfo_list.append(val)
        for key in handle.readline().strip().split("\t"):
            headers_list.append(key + "_2")
        for val in handle.readline().strip().split("\t"):
            bodyinfo_list.append(val)
        for key in handle.readline().strip().split("\t"):
            headers_list.append(key + "_single")
        for val in handle.readline().strip().split("\t"):
            bodyinfo_list.append(val)

    for key, val in zip(headers_list, bodyinfo_list):
        col_dict[key] = val
    return col_dict

def merge_statout(dir_name, flag, reads_type):
    '''merge each col_dict to a big dict'''
    path_list = []
    #os.chdir(dir_name)
    #for file in glob.glob("*.stat_out"):
    #    path_list.append(file)
    for file in os.listdir(dir_name):
        if file.endswith("stat_out"):
            path_list.append(os.path.join(dir_name, file))

    col_dict = {}
    statout_dict = {}
    statout_dict["headers"] = []
    statout_dict["bodyinfo"] = []

    for statout_path in path_list:
        statout_path = statout_path.strip()
        with open(statout_path, 'r') as handle:
            col_dict = parse_statout(handle, statout_path, flag, reads_type)
        statout_dict["bodyinfo"].append(col_dict)
    statout_dict["headers"] = col_dict.keys()
    sorted(statout_dict["headers"])

    #script_dir = os.path.abspath(sys.path[0])
    #os.chdir(script_dir)

    return statout_dict

def get_statout(dirname, statout_out, flag, reads_type):
    '''get stat_out'''
    statout_dict = merge_statout(dirname, flag, reads_type)
    write_csv(statout_out, statout_dict["headers"], statout_dict["bodyinfo"])

def main():
    '''main function'''
    parser = argparse.ArgumentParser(description='get fq stat info')
    parser.add_argument('--type', action='store', dest='readstype', help='reads type, SE or PE')
    parser.add_argument('--fqstatdir', action='store',
                        dest='fqstatdir', help='a dir contain all *fqStat.txt file')
    parser.add_argument('--fqstatout', action='store',
                        dest='fqstatout', help='output fqstat out info')
    parser.add_argument('--cleandir', action='store',
                        dest='cleandir', help="a dir contain all *clean.stat_out file")
    parser.add_argument('--cleanout', action='store',
                        dest='cleanout', help='output clean stat out info')
    parser.add_argument('--rmhostdir', action='store',
                        dest='rmhostdir', help="a dir contain all *rmhost.stat_out file")
    parser.add_argument('--rmhostout', action='store',
                        dest='rmhostout', help='output rmhost stat out info')
    args = parser.parse_args()

    get_fqstat(args.readstype, args.fqstatdir, args.fqstatout)
    get_statout(args.cleandir, args.cleanout, "clean", args.readstype)
    get_statout(args.rmhostdir, args.rmhostout, "rmhost", args.readstype)

if __name__ == "__main__":
    main()
