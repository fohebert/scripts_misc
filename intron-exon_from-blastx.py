#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                                  User Commands

\033[1mSYNOPSIS\033[0m
        Parses an output file from a blast on NCBI's website  (result file : txt)
        to  find  the exon-intron breaks. The sequences blasted against NCBI  can
        be of any type (e.g contigs from a de novo assembly) and they can be from
        gDNA or  cDNA.  The  script  verifies  the  positions on each contig that 
        blast for an NCBI sequence and writes it down in an output file. The pos-
        itions that blast for a protein in genbank (swissprot or nr) are necessa-
        rily coding regions and this is the premiss on which the program is based.

\033[1mUSAGE\033[0m
        %program <blast_output> <output>

\033[1mINPUT FILE - FORMAT\033[0m
        
        Blast output
        ************
        It has to be either the blast result file (.txt) from \033[1mNCBI's  website\033[0m  OR
        the blast result file from \033[1mBLAST+\033[0m (e.g local  blast,  swissprot,  nr,  nt
        or custom database) in FORMAT 0 (default).
        
        The program will look for the best blast result to output the exon-intron
        breaks. It looks for the best e-value of all the hits  for a given  query
        and uses its blasting positions for the final output.

\033[1mOUTPUT FILE - FORMAT\033[0m
        
        Contig          Exons
        contig_455      45-267; 897-1203
        contig_1434     670-830;
        ...             ...

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

queries = {}
query = ""
started = False

with open(blast_in, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        if line.startswith("Query="):
            
            if query != "":

                to_discard = set()
                if len(queries[query]) > 1:
                    
                    for hit in list(sorted(queries[query])):
                        
                        if hit not in to_discard:
                            current_evalue = queries[query][hit]["evalue"]
                            
                            for compared_hit in list(sorted(queries[query])):
                                if compared_hit != hit: 
                                
                                    if compared_hit not in to_discard:
                                        compared_evalue = queries[query][compared_hit]["evalue"]
                                        
                                        if current_evalue < compared_evalue:
                                            to_discard.add(compared_hit)
                                        elif current_evalue > compared_evalue:
                                            to_discard.add(hit)
                                        elif current_evalue == compared_evalue:
                                            if len(queries[query][hit]["positions"]) < len(queries[query][compared_hit]["positions"]):
                                                to_discard.add(compared_hit)
                                            elif len(queries[query][hit]["positions"]) > len(queries[query][compared_hit]["positions"]):
                                                to_discard.add(hit)
                                            else:
                                                to_discard.add(hit)

                    for hit_to_discard in to_discard:
                        queries[query].pop(hit_to_discard)
                
                elif len(queries[query]) == 1:
                    pass
            
            query = line.split("Query= ")[1]
            queries[query] = {}
        
        elif line.startswith("> "):
            
            hit = line.split("> ")[1]
            started = False
            queries[query][hit] = {}
        
        elif re.findall("Expect", line) != []:
            
            if started == False:
                evalue = float(line.split()[-1])
                queries[query][hit]["evalue"] = evalue
                queries[query][hit]["positions"] = set()
                started = True
            elif started:
                pass
        
        elif re.findall("Query  [0-9]*", line) != []:
            
            pos_sorting = []
            pos_sorting.append(int(line.split()[1]))
            pos_sorting.append(int(line.split()[-1]))
            
            pos1 = sorted(pos_sorting)[0]
            pos2 = sorted(pos_sorting)[1]
            
            for pos in range(pos1,pos2+1):
                queries[query][hit]["positions"].add(pos)

exon = set()
with open(out_file, "w") as out_f:
    
    out_f.write("Gene" + "\t" + "Exon positions")
    
    for query in list(sorted(queries)):
        out_f.write("\n" + query + "\t")
        
        if queries[query] != {}:
        
            for hit in queries[query]:
            
                lgth = len(queries[query][hit]["positions"])
                first_pos = 0
                for i in range(0,lgth):

                    if i == 0:
                        first_pos = sorted(queries[query][hit]["positions"])[i]
                    
                    elif sorted(queries[query][hit]["positions"])[i] == sorted(queries[query][hit]["positions"])[i-1] + 1:
                        pass
                    
                    elif sorted(queries[query][hit]["positions"])[i] > sorted(queries[query][hit]["positions"])[i-1] + 1:
                        out_f.write(str(first_pos) + "-" + str(sorted(queries[query][hit]["positions"])[i-1]) + ";")
                        first_pos = sorted(queries[query][hit]["positions"])[i]
                        
                    if i == (lgth - 1):
                        out_f.write(str(first_pos) + "-" + str(sorted(queries[query][hit]["positions"])[i]) + ";")
            
        elif queries[query] == {}:
            out_f.write("no exons found")
            
    print "\nFile processed\n"
