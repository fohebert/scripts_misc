#!/usr/bin/python

"""
                              User Commands

SYNOPSIS:
        Takes a blast+ output (format 6) from a blast of the contigs
        on the same contigs and creates an output file with the names
        of the contigs that blast against other ones. It returns a
        list of "redundant" contigs with the contigs to which they
        blast. 
        
        In other words, it verifies if there is contigs that are very 
        similar in a file of many apparently different contigs (from a 
        de novo assembly for example).

USAGE:
        %program <blast_fmt_6> <output_redundant_contigs>

CREDITS:
        Doctor Pants 2012 \m/

"""

import sys

try:
    blast_input = open(sys.argv[1], "rU")
    list_output = sys.argv[2]
except:
    print __doc__
    sys.exit(1)

contigs = []
redundant = {}

for line in blast_input:
    line = line.strip()
    if line != "":
        contig = line.split()[0]
        if contig in contigs:
            asso_contig = line.split()[1]
            if contig in redundant:
                redundant[contig].append(asso_contig)
            else:
                redundant[contig] = []
                redundant[contig].append(asso_contig)
        else:
            contigs.append(contig)

with open(list_output, "w") as out_f:
    for prim_contig in redundant:
        out_f.write(prim_contig)
        for sec_contig in redundant[prim_contig]:
            out_f.write("\t" + sec_contig)
        out_f.write("\n")
