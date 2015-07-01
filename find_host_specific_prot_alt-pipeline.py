#!/usr/bin/env python

"""
v.0.1                          User's Commands                           v.0.1

\033[1mDESCRIPTION\033[0m
    Takes the table output file from the script < blastp_parser_inclusive.py >
    and finds the contigs/proteins that  blasted  against  the  host  ONLY. In
    other words, it finds proteins that are exclusive to the host.

\033[1mUSAGE\033[0m
    %program <in_table> <out_file>

\033[1mCREDITS\033[0m
"""

import sys

try:
    in_table = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

host = raw_input("Name of the host as written in file: ")

count = 1
columns = {}
contigs = []
with open(in_table, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        if count == 1:
            num_col = len(line.split("\t"))
            col = 1
            while col <= num_col:
                columns[col] = line.split("\t")[col-1]
                col += 1
        
        elif count > 1:
            cont_name = line.split("\t")[0]
            col = 2
            h = 0
            others = 0
            while col <= len(line.split("\t")):
                if line.split("\t")[col-1] > 0:
                    if columns[col-1] == host:
                        h += 1
                    elif columns[col-1] != host:
                        others += 1    
                col += 1
            
            if h == 1 and others == 0:
                contigs.append(cont_name)
        
        count += 1

with open(out_file, "w") as out_f:
    for contig in contigs:
        out_f.write(contig + "\n")
        
print "\n\033[1mJob Done\033[0m\n"