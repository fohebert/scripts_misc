#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
                         \033[1mUser Commands\033[0m

\033[1mSYNOPSIS\033[0m:
        Extracts the sequences that are above a certain length.
        Requires BIO PYTHON to work.

\033[1mUASGE\033[0m:
        %program <fasta_input> <min_contig_lenght> <fasta_out>
        
\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys
from Bio import SeqIO

try:
    fasta_in = sys.argv[1]
    min_length = int(sys.argv[2])
    output = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

fasta_sequences = SeqIO.parse(open(fasta_in),"fasta")

with open(output, "w") as out_f:
    for seq in fasta_sequences:
        name = seq.id
        if len(seq) >= min_length:
            SeqIO.write([seq], out_f, 'fasta')
