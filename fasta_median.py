#!/usr/bin/env python

"""
\033[1mDESCRIPTION\033[0m
    Finds the median length of a given FASTA file.

\033[1mUSAGE\033[0m
    %program <fasta_in> <out_file>
    
    \033[1mNB\033[0m: <out_file> contains all sequence lengths
    
\033[1mCREDITS\033[0m
    Doc Pants 2014 \033[1m\m/\033[0m
"""

import sys
import numpy as np
from Bio import SeqIO

try:
    fasta_in = open(sys.argv[1], "rU")
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(1)

sequences = [(seq.id, seq.seq.tostring()) for seq in SeqIO.parse(fasta_in, "fasta")]

all_seq = []
for seq in sequences:
    all_seq.append(len(seq[1]))

med = np.median(all_seq)

print "\n\033[1mMedian length\033[0m: ", med, "\n"

with open(out_file, "w") as out_f:
    for seq in sorted(all_seq):
        out_f.write(str(seq) + "\n")