#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
v.0.1                  User's Commands               v.0.1

\033[1mDESCRIPTION\033[0m
    Takes  a FASTA  file  and generates overlapping k-mers (where "k"
    is specified in the command line arguments) with a sliding window
    of increment "x" (also specified by the user).

\033[1mUSAGE\033[0m
    %program <FASTA_in> <k-mer_length> <increment_length> <out_FASTA>

\033[1mCREDITS\033[0m
    DOCTOR PANTS \m/ 2013
"""

import sys
from Bio import SeqIO

try:
    in_fasta = open(sys.argv[1], "rU")
    k_mer = int(sys.argv[2])
    inc = int(sys.argv[3])
    out_file = sys.argv[4]
except:
    print __doc__
    sys.exit(0)

sequences = [(seq.id, seq.seq.tostring()) for seq in SeqIO.parse(in_fasta, "fasta")]

def k_mer_generator(seq, k_mer, inc): # s = sequence | k_mer = k-mer length | inc = increment length
    """Calls a function that can generate k-mers of specified length with specified increment
    based on a given a.a. sequence.
    """
    fragments = []
    s = seq[1]
    name = seq[0]
    num_frag = 1 # Count the number of k-mers per sequence to create their name
    count_pos_1 = 0
    count_pos_2 = k_mer
    
    # Loop that creates each k-mer for a given seq. Starts from pos 0 and creates a first fragment of
    # k-mer length. Then adds the increment length (given in the command line arguments) to the "counting
    # objects" (one for start and end position of the k-mer seq.) in order to generate fragments of
    # k-mer legnth with specified increment (ex. : 1 nucl.).
    while count_pos_2 <= len(s):
        frag = s[count_pos_1:count_pos_2]
        fragments.append([name + "_" + str("%04i" % num_frag), frag])
        count_pos_1 += 1
        count_pos_2 += 1
        num_frag += 1
    return fragments

count = 0
with open(out_file, "w") as out_f:
    for seq in sequences: # Goes through each sequence contained in the input FASTA file given
        fragments = k_mer_generator(seq, k_mer, inc) # For each seq., calls the 'k-mer generating function'
        for frag in fragments: # Elements in 'fragments' are tupples : (seq_name, k-mer_seq)
            # Writes down in output FASTA file the seq. of each k-mer generated for a given seq.
            if count == 0:
                out_f.write(">" + frag[0] + "\n" + frag[1])
            if count > 0:
                out_f.write("\n" + ">" + frag[0] + "\n" + frag[1])
            count += 1

print "\n\033[1mJob Done\033[0m\n"
