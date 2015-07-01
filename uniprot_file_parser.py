#!/usr/bin/env python

"""
v.1.0.                    User's Commands                    v.1.0

\033[1mDESCRIPTION\033[0m
    Takes the Flat Text output file from uniprot retrieve function
    on the website and parses the lines to extract all genes names
    found  in  the  file. Will create  an  input  file for GoMiner
    hightroughput online service.

\033[1mUSAGE\033[0m
    %program <input_swissprot> <output>
    
\033[1mCREDITS\033[0m
    Doc. Pants 2014 \m/
"""

import sys
import re

try:
    in_swissprot = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

with open(in_swissprot, "rU") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            line = line.strip()
            
            if line.startswith("GN"):
                if re.findall("Name=", line) != []:
                    gene = line.split("Name=")[1].split(";")[0]
                    
                    print gene
                    
                    if gene.find(",") != -1:
                        gene = gene.split(",")[0]
                        out_f.write(gene + "\n")
                    elif gene.find(",") == -1:
                        out_f.write(gene + "\n")
                elif re.findall("ORFNames=", line) != []:
                    gene = line.split("ORFNames=")[1].split(";")[0]
                    if gene.find(",") != -1:
                        gene = gene.split(",")[0]
                        out_f.write(gene + "\n")
                    elif gene.find(",") == -1:
                        out_f.write(gene + "\n")
                elif re.findall("OrderedLocusNames=", line) != []:
                    gene = line.split("OrderedLocusNames=")[1].split(";")[0]
                    if gene.find(",") != -1:
                        gene = gene.split(",")[0]
                        out_f.write(gene + "\n")
                    elif gene.find(",") == -1:
                        out_f.write(gene + "\n")

print "\033[1mJob Done\033[0m"