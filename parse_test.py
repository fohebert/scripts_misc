#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Prend un fichier txt avec des noms de séquences voulues et les cherche dans un fichier FASTA qui en
# contient une multitude et les ressort dans un fichier txt. Permet d'extraire rapidement un lot de 
# séquences (ex. : une centaine) d'un fichier qui en contient des milliers.

""" Takes a list of sequence names and extracts the sequences from a fasta file

    --> Usage : <input fasta file> <name file> <output file> <-- """

import sys
import re
from Bio import SeqIO

try:
    fasta_file = sys.argv[1]  # Input fasta file
    number_file = sys.argv[2] # Input interesting numbers file, one per line
    result_file = sys.argv[3] # Output fasta file

except:
    print __doc__
    sys.exit(0)

wanted = set()
with open(number_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
end = False
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")
    print "\n Extraction completed ! \n"
