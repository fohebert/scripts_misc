#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
\033[1mDescription\033[0m: 
    Takes a fasta file and eliminates any tabs in the sequence names

\033[1mUsage\033[0m:
    %program <FASTA_in> <FASTA_out>
"""

import sys
from Bio import SeqIO

prefix = raw_input("Enter gene prefix here : ")

try:
    fasta_in = open(sys.argv[1], "rU")
    fasta_out = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(fasta_in, "fasta"))

count = 0
with open(fasta_out, "w") as out_f:
    for seq in sequences:
        if count == 0:
            out_f.write(">" + prefix + seq[0] + "\n" + seq[1])
        else:
            out_f.write("\n" + ">" + prefix + seq[0] + "\n" + seq[1])
        
        count += 1

print "\nJob done\n"
