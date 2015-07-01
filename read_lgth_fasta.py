#!/usr/bin/env python

"""
                           User Commands

\033[1mSYNOPSIS\033[0m
        Takes a FASTA file and returns a list with the read lengths
        in a tab delimited text  file. One  read per line with ONLY
        the length for that read.

\033[1mUSAGE\033[0m
        %program <FASTA_file> <output>

\033[1mCREDITS\033[0m
        Doc The CrowPantartigan 2012 \m/
"""

import sys
from Bio import SeqIO

try:
    in_fasta = open(sys.argv[1], "rU")
    out = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(in_fasta, 'fasta'))

with open(out, "w") as out_f:
    for seq in sequences:
        out_f.write(str(len(seq[1])) + "\n")

print "\nDONE\n"
