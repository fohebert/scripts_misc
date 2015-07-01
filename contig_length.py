#!/usr/bin/python

"""

                       User Commands

DESCRIPTION
        Takes a FASTA file and returns in an output file
        the contig name and its length (tab text format).

USAGE
        %program <fasta_file> <result_file>
        
CREDITS
        Doctor Pants 2011 \m/

"""

import sys
from Bio import SeqIO

try:
    in_file = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

sequences = SeqIO.parse(open(in_file), "fasta")

with open(out_file, "w") as out_f:
    for seq in sequences:
        contig = seq.id
        length = str(len(seq))
        out_f.write(contig + "\t" + length + "\n")
    print "File processed !"
