#!/usr/bin/python

"""
                            User Commands

\033[1mSYNOPSIS\033[0m:

        Takes the output file from the script "bayescan_format.py"
        and generates a list of loci that must be discarded by the
        program bayescan, i.e the ones for which not all samples were
        successfully genotyped.

\033[1mUSAGE\033[0m
        
        %program <input_file> <output_file>

        input_file
        ******************************
        output from bayescan_format.py
        
        outpuf_file:
        ******************************
        file containing the loci to be discarded (one per line)
        
        \033[1m*IMPORTANT NOTE*\033[0m : the number of genes per pop.
        must be changed directly in the script.
        
\033[1mCREDITS\033[0m

        Doctor Pants 2011 \m/

"""

import sys
import re

try:
    in_file = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

# Specifies the maximum number of genes within one population
# (twice the amount of individuals in the pop if they are di-
# ploid)
max_num_genes = 6

loci_to_discard = set()

with open(in_file, "r") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            line = line.strip()
            if line != "":
                
                if line[0:6] == "[loci]":
                    continue
                if line[0:6] == "[popul":
                    continue
                if line[0:6] == "[pop]=":
                    continue
                if re.findall("^[1-9]", line):
                    countA = int(line.split()[3])
                    countB = int(line.split()[4])
                    
                    if countA + countB != max_num_genes:
                        loci_to_discard.add(line.split()[0] + "\n")
                        
        for l in list(sorted(loci_to_discard)):
            out_f.write(l)
                
        print "\nHow does it feel to be one of the beautiful people ?\n"
