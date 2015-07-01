#!/usr/bin/env python

"""
V.1.0                   User's Commands                v.1.0

\033[1mDESCRIPTION\033[0m
    Chooses the longest ORF for each sequence among multiple 
    sequences produced by  "getorf" (EMBOSS).
    
\033[1mUSAGE\033[0m
    %program <emboss_orf_file> <out_fasta_file>
        
\033[1mCREDITS\033[0m
    Doc Pants 2014
"""

import sys
from Bio import SeqIO
import re

try:
    emboss_in = open(sys.argv[1])
    out_fasta = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(emboss_in, "fasta"))

orfs = {}
kept = []
contigs = {}
for seq in sequences:
    
    # adds the specific ORF to the dict() containing all sequences
    # produced by EMBOSS getorf.  Will be  used later in  order to
    # retrieve the desired sequence, i.e. the longest one.
    contigs[seq[0]] = seq[1]
    
    # sets the original name of the contig (without the numbers
    # added to the sequence ID by EMBOSS getorf
    contig = seq[0].split("_")[0]
    orf_num = seq[0].split("_")[-1] # takes the orf number
    orf_lgth = len(seq[1]) # takes the orf length
    
    # adds the contig to the dict() orfs, which is used to store the
    # length info on each ORF found for each contig. ORFs are stored
    # based on the name of the contig from where they were extracted
    # and the number that was it assigned by EMBOSS, and finally,
    # their length.
    if contig not in orfs:
        
        orfs[contig] = {}
        orf_lgth = len(seq[1])
        orf_num = seq[0].split("_")[-1]
        orfs[contig][orf_num] = orf_lgth
        
    elif contig in orfs:
        
        orf_lgth = len(seq[1])
        orf_num = seq[0].split("_")[-1]
        orfs[contig][orf_num] = orf_lgth

to_keep = set()
with open(out_fasta, "w") as out_f:
    
    for contig in sorted(orfs):
        
        sorting = set()
        for orf_num in orfs[contig]:
            sorting.add(orfs[contig][orf_num])
        
        longest_orf = sorted(sorting)[-1]
        
        for orf_num in orfs[contig]:
            if orfs[contig][orf_num] == longest_orf:
                to_keep.add(contig + "_" + str(orf_num))
        
    for contig in list(sorted(contigs)):
    
        if contig in to_keep:
            out_f.write(">" + contig.split("_")[0] + "\n" + contigs[contig] + "\n")

print "\nDone\n"
