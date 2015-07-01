#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                           User Commands

\033[1mSYNOPSIS\033[0m
        This program parses the output file from BLAST+ in  order
        to retain only the sequences that have the best blast hit.
        The blast should be ORFs from multiple gene sequences bla-
        sted against a reference protein database (e.g translated
        transcriptome). The  program will only keep the ORFs that
        seem more  likely  to be exons and it will produce an out-
        put file that contains  the name  of the genes with their
        respective exons (their position).

\033[1mUSAGE\033[0m
        %program <blastplus_in> <output_file>

\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys
import re

try:
    blast_in = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

genes = {}
temporary = {}
current_gene = ""
queries = {}
hit_pos = set()
query_pos = set()
query_reverse = False
started = False
no_hits = False
end = False
hit = 1
with open(blast_in, "rU") as b_in:
    for line in b_in:
        line = line.strip()
        
        if line.startswith("Query="):
            
            if line.find("REVERSE") != -1:
                query_reverse = True
            elif line.find("REVERSE") == -1:
                query_reverse = False
            
            query_pos = set()
            gene = re.split("_[0-9]* \[", line.split("Query= ")[1])[0]
            query_pos.add(int(line.split("[")[1].split("]")[0].split()[0]))
            query_pos.add(int(line.split("[")[1].split("]")[0].split()[-1]))
            
            started = False
            end = False
            no_hits = False
            
            if current_gene == "":
                current_gene = gene
            
            elif current_gene == gene:
                pass
            
            elif current_gene != gene:
                
                to_discard = set()
                if temporary == {}:
                    genes[current_gene] = "No exons found"
                
                elif temporary != {}:
                    
                    genes[current_gene] = {}
                    
                    for current_hit in sorted(temporary):
                        
                        current_range = set()
                        if current_hit not in to_discard:
                            current_pos1 = sorted(temporary[current_hit]["pos"])[0]
                            current_pos2 = sorted(temporary[current_hit]["pos"])[-1]
                            for pos in range(current_pos1, current_pos2 + 1):
                                current_range.add(pos)
                            
                            for compared_hit in sorted(temporary):
                                if compared_hit == current_hit:
                                    pass
                                
                                elif compared_hit != current_hit and compared_hit not in to_discard:
                                    
                                    count = 0
                                    compared_pos1 = sorted(temporary[compared_hit]["pos"])[0]
                                    compared_pos2 = sorted(temporary[compared_hit]["pos"])[-1]
                                    for i in range(compared_pos1, compared_pos2 + 1):
                                        if i in current_range:
                                            count += 1
                                    
                                    if count >= 6:
                                        if temporary[current_hit]["evalue"] < temporary[compared_hit]["evalue"]:
                                            to_discard.add(compared_hit)
                                        elif temporary[current_hit]["evalue"] > temporary[compared_hit]["evalue"]:
                                            to_discard.add(current_hit)
                                        elif temporary[current_hit]["evalue"] == temporary[compared_hit]["evalue"]:
                                            if temporary[current_hit]["identities"] > temporary[compared_hit]["identities"]:
                                                to_discard.add(compared_hit)
                                            elif temporary[current_hit]["identities"] <= temporary[compared_hit]["identities"]:
                                                to_discard.add(current_hit)
                    
                    for hit_to_discard in sorted(to_discard):
                        temporary.pop(hit_to_discard)
                    
                    for hit in sorted(temporary):
                        exon_first_pos = sorted(temporary[hit]["pos"])[0]
                        exon = str(sorted(temporary[hit]["pos"])[0]) + "-" + str(sorted(temporary[hit]["pos"])[-1])
                        genes[current_gene][exon_first_pos] = exon
                    
                hit = 0
                current_gene = gene
                temporary = {}
            
            hit += 1
            
        elif line.find("No hits found") != -1:
            no_hits = True
        
        elif line.startswith("> ") and started:
            if end == False:
                end = True
            elif end:
                pass
        
        elif line.startswith("Score =") and started == False:
            
            started = True
            temporary[hit] = {}
            temporary[hit]["evalue"] = 0
            temporary[hit]["evalue"] += float(line.split("Expect")[-1].split()[-1])
            temporary[hit]["identities"] = 0
            temporary[hit]["pos"] = set()
            
        elif line.startswith("Identities") and started and end == False:
            hit_pos = set()
            temporary[hit]["identities"] += int(line.split("(")[1].split(")")[0].split("%")[0])
        
        elif re.findall("Query  [0-9]*", line) != [] and started and end == False:
            
            pos1 = int(line.split()[1])
            pos2 = int(line.split()[-1])
            
            hit_pos.add(pos1)
            hit_pos.add(pos2)
            
        elif line.startswith("Effective search") and no_hits == False:
            
            query_pos1 = sorted(query_pos)[0]
            query_pos2 = sorted(query_pos)[-1]
            hit_pos1 = sorted(hit_pos)[0]
            hit_pos2 = sorted(hit_pos)[-1]
            
            if query_reverse == False:
                temporary[hit]["pos"].add(query_pos1 + (hit_pos1 - 1))
                temporary[hit]["pos"].add(query_pos2 - ((query_pos2 - query_pos1 + 1) - hit_pos2))
                
            if query_reverse:
                temporary[hit]["pos"].add(((query_pos2 - query_pos1 + 1) - hit_pos2) + query_pos1)
                temporary[hit]["pos"].add(query_pos2 - (hit_pos1 - 1))
        
with open(out_file, "w") as out_f:
    
    out_f.write("Gene" + "\t" + "Exon positions")
    
    for gene in list(sorted(genes)):
        out_f.write("\n" + gene + "\t")
        
        if genes[gene] == "No exons found":
            out_f.write("No exons found")
        
        elif genes[gene] != "No exons found":
        
            for exon_first_pos in sorted(genes[gene]):
                out_f.write(genes[gene][exon_first_pos] + ";")
            
print "\nFile processed and output file successfully created\n"
