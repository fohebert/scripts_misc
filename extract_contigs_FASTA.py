#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
                         \033[1mUser Commands\033[0m

\033[1mSYNOPSIS\033[0m:
        Extracts the sequences of contigs (FASTA file) that are listed in
        a text file, one name per line, and outputs them in a FASTA file 
        (output file).

\033[1mUASGE\033[0m:
        %program <fasta_input> <wanted_names> <fasta_out>
        
\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys
from Bio import SeqIO

try:
    fasta_in = open(sys.argv[1], "rU")
    wanted_contigs = sys.argv[2]
    fasta_out = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

wanted = set()
with open(wanted_contigs, "r") as in_f:
    for line in in_f:
        line = line.strip()
        if line != "":
            wanted.add(line)

sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(fasta_in, "fasta"))

with open(fasta_out, "w") as out_f:
    for seq in sequences:
        name = seq[0]
        if name in wanted:
            out_f.write(">" + name + "\n" + seq[1] + "\n")

contigs_new_fasta = set()
with open(fasta_out, "rU") as final_f:
        
    for line in final_f:
        line = line.strip()
        if line.startswith(">"):
            contigs_new_fasta.add(line.split(">")[1])
    
    for name in wanted:
        if name not in contigs_new_fasta:
            print "Missing :", name
    
    if len(contigs_new_fasta) == len(wanted):
        print "\nAll contigs have been found and successfully extracted !\n"
    elif len(contigs_new_fasta) != len(wanted):
        print "Problem with the names in <reference_fasta> and <wanted_text> file\nExtraction incomplete"
