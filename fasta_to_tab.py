#!/usr/bin/python

from Bio import SeqIO
import sys

fasta_file = sys.argv[1]
result_file = sys.argv[2]

sequences = SeqIO.parse(open(fasta_file), 'fasta')

with open(result_file, "w") as f:
    for seq in sequences:
        SeqIO.write([seq], f, 'tab')
    print "Done !"
