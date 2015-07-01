#! /usr/bin/env python
# -*- coding : utf-8 -*-

"""
v.0.1                 User's Commands              v.0.1

\033[1mDESCRIPTION\033[0m
    Extracts sequences from a FASTA file that  are consi-
    dered  as  "mimicry sequences" based on BLAST results 
    of these sequences against one  or several  reference
    proteomes (e.g. proteome from a  parasite host). If a 
    given  query  shows  high similarity with any control 
    sequence, it is written in a FASTA output file.
    
    Also produces a table output file with information on
    the identity score of each query sequence.
    
    \033[1mNB : BLASTp result file has to be in fmt 0.\033[0m

\033[1mUSAGE\033[0m
    %program <BLAST_results> <FASTA_in> <FASTA_out> <table>

\033[1mCREDITS\033[0m
    DOC PANTS \m/ 2013
"""

# Importing modules
import sys
from Bio import SeqIO
import re

try:
    in_blast = sys.argv[1]
    in_fasta = sys.argv[2]
    out_fasta = sys.argv[3]
    out_table = sys.argv[4]
except:
    print __doc__
    sys.exit(2)

# BLASTp identity thresholds placed in a dict() as a reference for later comparison
# Thresholds empirically pre-defined by Ludin et al. 2011
t_denom = (10,11,12,13,14)
t_num = (10,10,11,12,12)
thresholds = dict(zip(t_denom,t_num))

k_mers = {} # global dict() containing the identity scores of all k_mers with hits
start = False
new_hit = False
query = ""
to_add = set()
with open(in_blast, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        # When there is a new query, retains its name and tells the program
        # that the investigation of a new query has started.
        if line.startswith("Query=") and start == False:
            query = line.split("Query= ")[-1]
            start = True
        
        # If it turns out that the current query has no hits found, stops everything
        # until the next query.
        elif re.findall("No hits found", line) != [] and start:
            query = ""
            start = False
        
        # If the current query has at least one hit, it tells the program that
        # it can starts looking at identity scores. If this is a new query not
        # contained in the global dict(), it adds the query to the global dict
        elif line.startswith("> ") and start:
            if query not in k_mers:
                k_mers[query] = []
                new_hit = True
            elif query in k_mers:
                new_hit = True
        
        # When it encounters the first score (best one) of a given hit, appends
        # it to the global dict() in the corresponding query "key". Will do this
        # for each best score of every hit for a given query.
        elif line.startswith("Identities") and start and new_hit:
            ident_num = line.split()[2].split("/")[0]
            ident_denom = line.split()[2].split("/")[-1]
            k_mers[query].append((ident_num, ident_denom))
            new_hit = False
        
        # When it has finished, it stops everything and it's ready to start all
        # over again for a new query.
        elif line.startswith("Effective search space") and start:
            query = ""
            start = False

with open(out_fasta, "w") as out_f:
    with open(out_table, "w") as out_t:
        # Looks in the global dict() and verifies, for each query, if it has at least one
        # identity score above pre-defined thresholds. If there is at least one identity
        # score above one of the thresholds, appends the query name to the set() containing
        # names to be discarded.
        for k_mer in sorted(k_mers):
            above_threshold = False
            for identity_score in k_mers[k_mer]:
                if int(identity_score[1]) >= 10:
                    if above_threshold == False and int(identity_score[0]) >= thresholds[int(identity_score[1])]:
                        to_add.add(k_mer)
                        out_t.write(k_mer + "\t" + identity_score[0] + "\t" + identity_score[1] + "\n")
                        above_threshold = True
        
        print "\n\033[1m'to_add' length: \033[0m", len(to_add)
        
        with open(in_fasta, "rU") as in_fas:
            # Generates a list object ("sequenceS") containing tuples with all sequences with
            # their respective name (ID)
            sequences = [(seq.id, seq.seq.tostring()) for seq in SeqIO.parse(in_fas, "fasta")]
        
            # Writes in the FASTA output file the sequences that passed the thresholds.
            count = 0
            for seq in sequences:
                if seq[0] in to_add:
                    if count == 0:
                        out_f.write(">" + seq[0] + "\n" + seq[1])
                    elif count > 0:
                        out_f.write("\n" + ">" + seq[0] + "\n" + seq[1])
                count += 1

print "\n\033[1mJob Done\033[0m\n"
