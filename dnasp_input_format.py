#!/usr/bin/env python

"""
DESCRIPTION
        Takes a FASTA file and formats it in a dnaSP input file.

USAGE
        %program <FASTA_in> <dnasp_input_file>
"""

import sys
import re

try:
    fasta_in = sys.argv[1]
    dnasp_in = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

seq_lgth = 0
lgth_taken = False
second_a = False
seq1 = ""
seq2 = ""
with open(fasta_in, "rU") as in_f:
    with open(dnasp_in, "w") as out_f:
        for line in in_f:
            line = line.strip()
            
            if line.startswith(">"):
            
                if line.find("_A1_") != -1:
                    seq_lgth = 0
                    lgth_taken = False
                    second_a = False
                
                elif line.find("_A2_") != -1:
                    second_a = True
                
            elif line.startswith(">") == False:
                    
                if lgth_taken == False:
                    seq1 = ""
                    lgth_taken = True
                    seq_lgth = len(line)
                    seq1 += line
                    
                elif lgth_taken and second_a:
                    seq2 = ""
                    seq2 += line
                    
                    out_f.write("2" + "\t" + str(seq_lgth) + "\n" + "\n" + "allele1" + "\n" + "\n" + seq1)
                    out_f.write("\n" + "\n" + "allele2" + "\n" + "\n" + seq2 + "\n" + "\n")
