#!/usr/bin/python

"""
v.0.1                      User's Commands                     v.0.1

\033[1mDESCRIPTION\033[0m
    Takes a concatenated blast result file (multiple blast results of
    format "0" concatenated) and returns in a FASTA output file  the
    sequences with an e-value > user-determined threshold.

\033[1mUSAGE\033[0m
    %program <blast_results> <FASTA_all_seq> <output_FASTA_filtered>

\033[1mCREDITS\033[0m
    Doc PANTS \m/ 2013
"""

import sys
from Bio import SeqIO
import re

try:
    in_blast = sys.argv[1]
    in_fasta = open(sys.argv[2], "rU")
    out_fasta = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

threshold = float(raw_input("e-value threshold = "))

sequences = [(seq.id, seq.seq.tostring()) for seq in SeqIO.parse(in_fasta, "fasta")]

proteins = {} # Dict containing the info on all proteins, i.e. protein name + e-values
to_discard = set() # proteins to be discarded because they have e-value <= 1e-10
start = False
score = False
protein = ""
with open(in_blast, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        # When a new query is found, protein name is retained
        if line.startswith("Query=") and start == False:
            protein = line.split("Query= ")[-1]
            start = True # indicates that the process of parsing a new query has begun
        
        # If the query we are currently looking at has no hit, it resets everything to the begining
        # i.e. empties the object containing the protein name and returns to start = False
        elif re.findall("No significant hits", line) != [] and start:
            protein = ""
            start = False
        
        # If the current query has some hits, the protein name is added in the dict.
        elif line.startswith("Sequences producing significant alignments:") and start:
            score = True
            
            if protein not in sorted(proteins):
                proteins[protein] = set() # creates a set() that will contain e-values for that protein
        
        # Takes the e-value for the best hit (i.e. the first in the list, see output file to confirm)
        # and puts it in the set() created for that protein (query). Then, it resets everything: back
        # to the start with protein name empty and score/start = False.
        elif re.findall("[a-z]_[a-z]*_[a-z]*", line) != [] and start and score:
            e_value = float(line.split()[-1])
            proteins[protein].add(e_value)
            score = False
            start = False
            protein = ""
    
    # Goes through the dict() and sorts all e-values (smaller -> larger) found for each protein. If
    # the smaller e-value (i.e. the "best" e-value) <= 1e-10 (threshold for conserved proteins), it
    # adds the protein name to the object containing the names to discard.
    for protein in list(sorted(proteins)):
        if sorted(proteins[protein])[0] <= threshold:
            to_discard.add(protein)
    
    print "Number of sequences to discard: ", len(to_discard)

# Goes through the FASTA file given in arguments in command line and if the name of the protein is
# in the object "to_discard", it does not write it in the output file, otherwise, the sequence is
# added to the output file, which contains only the sequences that had an e-value > 1e-10.
count = 0
with open(out_fasta, "w") as out_f:
    for seq in sequences:
        if seq[0] not in sorted(to_discard):
            if count == 0:
                out_f.write(">" + seq[0] + "\n" + seq[1])
            if count > 0:
                out_f.write("\n" + ">" + seq[0] + "\n" + seq[1])
        count += 1

print "\n\033[1mJob done\033[0m\n"
