#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
\033[1mUSAGE\033[0m
        %program <hetero_file> <list_outliers> <output_all> <output_ouliers>
"""

import sys

try:
    in_hetero = sys.argv[1]
    list_outliers = sys.argv[2]
    output_all = sys.argv[3]
    output_outliers = sys.argv[4]
except:
    print __doc__
    sys.exit(2)

outliers = {}
with open(list_outliers, "rU") as list_o:
    for line in list_o:
        line = line.strip()
        
        gene = line.split("\t")[0]
        pos = line.split("\t")[1]
        
        if gene not in outliers:
            outliers[gene] = set()
            outliers[gene].add(pos)
        
        elif gene in outliers:
            outliers[gene].add(pos)

with open(in_hetero, "rU") as in_h:
    with open(output_all, "w") as out_all:
        with open(output_outliers, "w") as out_outliers:
            for line in in_h:
                
                if line.startswith("Gene") == False:
                    
                    gene = line.split("\t")[0]
                    pos = line.split("\t")[1]
                    Ho_d = line.split("\t")[2]
                    Ho_n = line.split("\t")[3]
                    
                    if gene in list(sorted(outliers)):
                        if pos in list(sorted(outliers[gene])):
                            out_outliers.write(gene + "\t" + pos + "\t" + Ho_d + "\t" + Ho_n)
                        
                        elif pos not in list(sorted(outliers[gene])):
                            out_all.write(gene + "\t" + pos + "\t" + Ho_d + "\t" + Ho_n)
                    
                    elif gene not in list(sorted(outliers)):
                        out_all.write(gene + "\t" + pos + "\t" + Ho_d + "\t" + Ho_n)
