#!/usr/bin/python

"""
SYNOPSIS:

    This program takes CLC's SNP detection output file
    and extracts only the SNPs from wanted contigs
    among the result file. The names of the wanted
    contigs are given in a text file, one name per line.

USAGE:

    %program <CLC SNP table> <contig names> <output>
    
CREDITS:
    
    Doc Pants 2011 \m/
"""

import sys

try:
    SNP_detection = sys.argv[1]
    contig_names = sys.argv[2]
    output = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

wanted = set()

with open(contig_names) as f:
    for line in f:
        l = line.strip()
        if l != "":
            wanted.add(l)

with open(SNP_detection) as in_f:
    with open(output, "w") as out_f:
        for line in in_f:
            contig = line.split("\t")[0]
            if contig in wanted:
                out_f.write(line)
        print "Done !"
