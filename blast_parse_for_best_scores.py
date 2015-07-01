#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                              User's Commands

\033[1mDESCRIPTION\033[0m
    Takes a BLAST result file (format 0) and extracts the single best score
    of each query (if it has at least one hit and in that case, it extracts
    the only hit available). It produces an output file with  the gene name
    in one column and in the other column,  it gives the  log2(best_score).
    This can be used to produce a  scatter plot  graph with the log2 of the
    best scores for a given proteome blasted on another proteome (see Ludin
    et al. 2011).

\033[1mUSAGE\033[0m
    %program <input_blast_result> <output_table>

\033[1mCREDITS\033[0m
    Doctor Pant 2013 \m/
"""

import sys
import math
import re

try:
    in_blast_result = sys.argv[1]
    out_table = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

count = 0
score = False
start = False
protein_name = ""
with open(in_blast_result, "rU") as in_f:
    with open(out_table, "w") as out_f:
        for line in in_f:
            line = line.strip()
            
            if line.startswith("Query=") and start == False:
                protein_name = line.split("Query= ")[-1]
                start = True
            
            elif re.findall("No hits found", line) != [] and start:
                start = False
                protein_name = ""
            
            elif re.findall("(Bits)", line) != [] and start:
                score = True
                if count == 0:
                    out_f.write(protein_name)
                if count > 0:
                    out_f.write("\n" + protein_name)
                count += 1
            
            elif re.findall("^Score =", line) != [] and start and score:
                score = float(line.split()[2].split(" bits")[0])
                out_f.write("\t" + str(math.log(score,2)))
                start = False
                score = False

print "\n\033[1mJob done\033[0m\n"