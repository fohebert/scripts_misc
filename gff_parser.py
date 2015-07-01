#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
\033[1mSYNOPSIS\033[0m
        Takes a .gff file from the program AUGUSTUS (exon prediction)
        and extracts the information concerning the exon positions of
        the  sequences  processed through AUGUSTUS. Creates an output
        file  containing  the name of the genes or the sequences pro-
        cessed and the  positions  on these sequences that correspond
        to putative exons.

\033[1mUSAGE\033[0m
        %program <GFF_file> <output_file>
        
\033[1mOUTPUT FILE - FORMAT\033[0m

                _________file starts below the line____________
                
                Gene                Exon Positions
                IGF-B               10-452;982-1223
                HSP-70              392-662;782-1912;2120-2532
                ...                 ...
                
                ___________file ends above the line____________

\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys
import re

try:
    in_gff = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

start = False
no_exons = False
exons_found = 0
exon_lgth = []
with open(in_gff, "rU") as i_gff:
    with open(out_file, "w") as out_f:
        
        out_f.write("Gene" + "\t" + "Exon positions")
        
        for line in i_gff:
            line = line.strip()
            
            if line.startswith("# ----- prediction on sequence number"):
                start = False
                gene = line.split("name = ")[-1].split(") -----")[0]
                out_f.write("\n" + gene + "\t")
            
            elif line.startswith("# Predicted genes for sequence number") and start == False:
                start = True
            
            elif line.find("\tCDS") != -1 and start:
                exon_pos1 = line.split()[3]
                exon_pos2 = line.split()[4]
                lgth = int(exon_pos2) - int(exon_pos1) + 1
                exon_lgth.append(lgth)
                out_f.write(exon_pos1 + "-" + exon_pos2 + ";")
                exons_found += 1
                
            elif line.startswith("# (none)") and start:
                out_f.write("\t" + "No exons found")
            
            elif line.startswith("# protein sequence") and start:
                start = False

mean_exon_lgth = float(sum(exon_lgth))/float(len(exon_lgth))
print "\n\033[1mMean exon length =\033[0m ", mean_exon_lgth, "(",sorted(exon_lgth)[0],"-",sorted(exon_lgth)[-1],")"
print "\033[1mNumber of exons found\033[0m: ", exons_found
print "\nDone\n"
