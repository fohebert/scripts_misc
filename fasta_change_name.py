#!/usr/bin/env python

"""
\033[1mDESCRIPTION\033[0m
    Takes a FASTA file and changes the name when there is
    multiple words with spaces. Keeps only the first word
    and  discards  everything  else that  is separated by
    either a tab or a space.

\033[1mUSAGE\033[0m
    %program <input_fasta> <output_fasta>
"""

import sys
from Bio import SeqIO

try:
    in_fasta = open(sys.argv[1], "rU")
    out_fasta = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

sequences = [(seq.id, seq.seq.tostring()) for seq in SeqIO.parse(in_fasta, "fasta")]

count = 0
with open(out_fasta, "w") as out_f:
    for seq in sequences:
        new_name = seq[0].split()[0]
        if count == 0:
            out_f.write(">" + new_name + "\n" + seq[1])
        else:
            out_f.write("\n" + ">" + new_name + "\n" + seq[1])
        count += 1

print "\n\033[1mJob Done\033[0m\n"
