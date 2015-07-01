#! /usr/bin/python

"""
                    User's Commands

\033[1mDESCRIPTION\033[0m
    Concatenates score columns from 2 different files in
    a single file that will be used to produce R graphs.

\033[1mUSAGE\033[0m
    %program <scores_host> <score_comparison> <out_file>

\033[1mCREDITS\033[0m
    DOC FUCKING PANTS \m/ 2013
"""

import sys

try:
    in_host = sys.argv[1]
    in_comparison = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(2)

scores = {}
both_present = set()
with open(in_host, "rU") as in_f1:
    with open(in_comparison, "rU") as in_f2:
        for line in in_f1:
            
            protein = line.split()[0]
            scores[protein] = {}
            line = line.strip()
            score1 = line.split()[-1]
            scores[protein]["species1"] = score1

        print len(scores)
        
        for line in in_f2:
            
            protein = line.split()[0]
            if protein in sorted(scores):
                both_present.add(protein)
                line = line.strip()
                score2 = line.split()[-1]
                scores[protein]["species2"] = score2
    
    with open(out_file, "w") as out_f:
        start = 0
        for protein in sorted(scores):
            if protein in both_present:
                if start == 0:
                    out_f.write(protein + "\t" + scores[protein]["species1"] + "\t" + scores[protein]["species2"])
                elif start > 0:
                    out_f.write("\n" + protein + "\t" + scores[protein]["species1"] + "\t" + scores[protein]["species2"])
                start += 1

print "\n\033[1mJob done\n\033[0m"
