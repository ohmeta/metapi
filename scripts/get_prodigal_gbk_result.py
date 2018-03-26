#!/usr/bin/python
# email: zhujie@genomics.cn
# license: GPL V3
import re

gbklist = "./gene.coordinate.gbk.pathlist.new"
out = open("./gene.coordinate.stat.out.new", 'w')
out.write("ID\tpartial=00\tpartial=01\tpartial=10\tpartial=11\ttotal_len\ttotal_num\tavg_length\n")

with open(gbklist, 'r') as path_handler:
    for gbkpath in path_handler:
        genenum = {}
        gene_total_len = 0
        gene_total_num = 0
        gene_avg_len = 0
        partial = ['partial=00', 'partial=01', 'partial=10', 'partial=11']
        genenum['partial=00'] = 0
        genenum['partial=01'] = 0
        genenum['partial=10'] = 0
        genenum['partial=11'] = 0

        with open(gbkpath.strip(), 'r') as gbk_handler:
            first = next(gbk_handler)
            id = re.search(r'(.*?)seqhdr="(CL\d+_L\d+_\d+)_scaffold(.*)', first).group(2)
            genenum['id'] = id
            for line in gbk_handler:
                for tag in partial:
                    if re.search(tag, line):
                        genenum[tag] += 1
                        gene_total_num += 1
                if re.search("CDS\s+(complement\()?(<)?(\d+)\.\.(>)?(\d+)(\))", line):
                    len = re.search("CDS\s+(complement\()?(<)?(\d+)\.\.(>)?(\d+)(\))", line)
                    gene_total_len += int(len.group(5)) - int(len.group(3)) + 1
            
            gene_avg_len = round(float(gene_total_len) / float(gene_total_num), 6)

            out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n" % (
                    genenum['id'],
                    genenum['partial=00'],
                    genenum['partial=01'],
                    genenum['partial=10'],
                    genenum['partial=11'],
                    gene_total_len,
                    gene_total_num,
                    gene_avg_len))

out.close()
