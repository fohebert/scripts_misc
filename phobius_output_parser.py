#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
v.0.1              User's Commands               v.0.1

\033[1mDESCRIPTION\033[0m
    Phobius output file parser: goes through  the file
    and extracts the position of the N-terminal signal
    peptide  of each sequence. Using this information,
    it extracts the a.a. sequence in the corresponding
    FASTA file \033[1mWITHOUT\033[0m the signal.

\033[1mUSAGE\033[0m
    %program <phobius_output_file> <FASTA> <output>

\033[1mCREDITS\033[0m
    DOCTOR PANTS \m/ 2013
"""

import sys
from Bio import SeqIO
import re

try:
    in_phobius = sys.argv[1]
    in_fasta = open(sys.argv[2], "rU")
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

sequences = [(seq.id, seq.seq.tostring()) for seq in SeqIO.parse(in_fasta, "fasta")]

start = False
signals = {}
protein = ""
with open(in_phobius, "rU") as in_p:
    for line in in_p:
        line = line.strip()
        
        # This is wehere the analysis of one sequence begings, by turning the object "start" on
        if line.startswith("ID"):
            protein = line.split()[-1]
            start = True
        
        # If it has started and if it finds the word "SIGNAL", it means there is an N-terminal
        # signal in that particular sequence. It enters the position of the end of that signal
        # in the dict() named "signals" (key = protein name, value = end of N-termina signal).
        elif line.find("SIGNAL") != -1 and start:
            signal_end = int(line.split()[-1])
            signals[protein] = signal_end
            if int(line.split()[2]) > 5:
                print "Signal not at begining: ", protein, line.split()[2], line.split()[-1]
        
            start = False
            protein = ""

count = 0
with open(out_file, "w") as out_f:
    for seq in sequences:
        
        # Goes through all the sequences from the original FASTA (i.e. without discarded
        # signals) and if the sequence is in the dict() "signals" (whicih means it has
        # a signal to be discared from the sequence) it writes down in the output FASTA
        # file the name of the sequence and the corresponding a.a. sequence, without the
        # the signal: starts the sequence from the end of the N-terminal signal.
        if seq[0] in sorted(signals):
            new_seq = seq[1][signals[seq[0]]:]
            if count == 0:
                out_f.write(">" + seq[0] + "\n" + new_seq)
            if count > 0:
                out_f.write("\n" + ">" + seq[0] + "\n" + new_seq)
        
        if seq[0] not in sorted(signals):
            if count == 0:
                out_f.write(">" + seq[0] + "\n" + seq[1])
            if count > 0:
                out_f.write("\n" + ">" + seq[0] + "\n" + seq[1])
        
        count += 1

print "\n\033[1mJob Done\033[0m\n"
