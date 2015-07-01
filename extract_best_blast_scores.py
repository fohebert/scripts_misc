#!/usr/bin/env python

"""
v.1.0                            User's Commands                            v.1.0

\033[1mDESCRIPTION\033[0m
    Parses 2 blast output files (fmt0) to keep only the best score for each query.
    Queries in blast output files were  blasted  against  the  same database, e.g.
    two different species (parasite for instance),  each  blasted against the same 
    species (host for instance). The program will produce pairs of best scores for
    the same target sequence. If a given query in one species blas ted against the
    exact same target sequence  than another query in the  second  species, then a
    pair with their respetive best score will be produced. Each pair produced that
    way  will be  a  single point on scatter plot in R representing the comparison
    of the best blast scores between two species when blasted on the same species.
    
\033[1mUSAGE\033[0m
    %program <blast_results_1> <blast_results_2> <output_file>
    
\033[1mCREDITS\033[0m
    Doc Pants 2014 \m/
"""

import sys
import re
import math

try:
    in_blast1 = sys.argv[1]
    in_blast2 = sys.argv[2]
    output_file = sys.argv[3]
except:
    print __doc__
    sys.exit(2)

scores = {}
start = False
hit = False
hit_id = ""

# Will get, for each query, the first hit with the total bit score for that hit
# Will then store the hid id in the dict() named "score", with its respective
# score as "the first score". The common aspect between both blast output files
# is the database on which the queries have been blasted, so hit IDs will serve
# as keys for the "scores" dict().
with open(in_blast1, "rU") as in_f1:
    for line in in_f1:
        line = line.strip()
        
        if line.startswith("Query=") and start == False:
            start = True
        
        elif re.findall("No hits found", line) != [] and start:
            start = False
        
        elif line.startswith(">") and start and hit == False:
            hit_id = line.split("> ")[-1]
            if hit_id not in scores:
                scores[hit_id] = {}
                scores[hit_id]["score1"] = []
            hit = True
        
        elif re.findall("bits \(", line) != [] and start and hit:
            score1 = float(line.split("Score = ")[1].split(" bits")[0])
            scores[hit_id]["score1"].append(str(math.log(score1,2)))
            hit = False
            start = False
            hit_id = ""

start = False
hit = False
hit_id = ""
# Does the same thing as above, but with the second blast output file.
# Will add the second score (named "score2") to the corresponding hit ID
# in the dict() "scores".
with open(in_blast2, "rU") as in_f2:
    for line in in_f2:
        line = line.strip()
        
        if line.startswith("Query=") and start == False:
            start = True
        
        elif re.findall("No hits found", line) != [] and start:
            start = False
        
        elif line.startswith("> ") and start and hit == False:
            hit_id = line.split("> ")[-1]
            hit = True
        
        elif re.findall("bits \(", line) != [] and start and hit:
            score2 = float(line.split("Score = ")[1].split(" bits")[0])
            
            if hit_id in scores:
                if len(scores[hit_id]) == 1:
                    scores[hit_id]["score2"] = []
                    scores[hit_id]["score2"].append(str(math.log(score2,2)))
                elif len(scores[hit_id]) == 2:
                    scores[hit_id]["score2"].append(str(math.log(score2,2)))
        
        elif line.startswith("Effective search space used:") and start and hit:
            start = False
            hit = False
            hit_id = ""

count = 0
with open(output_file, "w") as out_f:
    for hit_id in sorted(list(scores)):
        
        if len(scores[hit_id]) == 2:
            
            if len(scores[hit_id]["score1"]) == len(scores[hit_id]["score2"]):
                
                num_scores = len(scores[hit_id]["score1"]) - 1
                while num_scores >= 0:
                    if count == 0:
                        out_f.write(hit_id + "\t" + scores[hit_id]["score1"][num_scores] + "\t" + scores[hit_id]["score2"][num_scores])
                    elif count > 0:
                        out_f.write("\n" + hit_id + "\t" + scores[hit_id]["score1"][num_scores] + "\t" + scores[hit_id]["score2"][num_scores])
                    count += 1
                    num_scores -= 1
                
            elif len(scores[hit_id]["score1"]) > len(scores[hit_id]["score2"]):
                
                num_scores = len(scores[hit_id]["score2"]) - 1
                while num_scores >= 0:
                    if count == 0:
                        out_f.write(hit_id + "\t" + scores[hit_id]["score1"][num_scores] + "\t" + scores[hit_id]["score2"][num_scores])
                    elif count > 0:
                        out_f.write("\n" + hit_id + "\t" + scores[hit_id]["score1"][num_scores] + "\t" + scores[hit_id]["score2"][num_scores])
                    count += 1
                    num_scores -= 1
                
            elif len(scores[hit_id]["score1"]) < len(scores[hit_id]["score2"]):
                
                num_scores = len(scores[hit_id]["score1"]) - 1
                while num_scores >= 0:
                    if count == 0:
                        out_f.write(hit_id + "\t" + scores[hit_id]["score1"][num_scores] + "\t" + scores[hit_id]["score2"][num_scores])
                    elif count > 0:
                        out_f.write("\n" + hit_id + "\t" + scores[hit_id]["score1"][num_scores] + "\t" + scores[hit_id]["score2"][num_scores])
                    count += 1
                    num_scores -= 1

print "\n\033[1mJob Done\033[0m\n"