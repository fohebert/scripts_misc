#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

SYNOPSIS:

    Takes the output of Blast2Go ("annot_GOs_<number>_<number>.txt")
    and creates an output file containing one contig per line with
    the corresponding GO information found by Blast2Go (tab separa-
    ted).
    
    Output format:
    
    Contig   Description   Molecular function   Biological Process   Cellular component
    
    NB : the input file has to be parsed manually in order to keep
         only 1 result per GO category found for each contig.

USAGE:

    %program <Blast2Go output> <contig_description> <out_file>
    
    contig_description format:
    
    Blast ID for each contig, one contig ID per line.
    
    \033[1mWARNING !!! DOESN'T WRITE THE INFO FOR THE LAST CONTIG, MUST
    MODIFY THIS DETAIL !\033[0m

CREDITS:

    Doctor Pants 2011 \m/

"""

import sys
import re

try:
    raw_format = sys.argv[1]
    contig_description = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit()

# Generates a set containing all the contig descriptions (blast results from Blast2Go).
description = set()
with open(contig_description) as f:
    for line in f:
        l = line.strip()
        if l != "":
            description.add(l)


contig = ""
treatment = ""
blast = ""
function = {}
process = {}
compartment = {}
count = 0

with open(raw_format) as in_file:
    with open(out_file, "w") as out_file:
        
        # Writes out the first line of the output file (column titles)
        out_file.write("Contig" + "\t" + "Description" + "\t" + "Molecular function" + "\t" + 
        "Biological process" + "\t" + "Cellular component" + "\n")
        
        for line in in_file:
        
            count += 1
            line = line.strip()

            # Writes the GO results for each contig after the dictionaries have been filled
            if line.split("\t")[1] in description and count > 1:
                
                # Fills in the empty dictinaries
                if function == {}:
                    function = {"unkown":"-"}
                if process == {}:
                    process = {"unkown":"-"}
                if compartment == {}:
                    compartment = {"unkown":"-"}
                
                # Writes in the output file the GO results, one dictionary at a time (one
                # dictionary per GO category)
                category = ""
                number = ""
                for k, v in function.iteritems():
                    category = k
                    number = v
                    out_file.write(contig + "\t" + blast + "\t" + category + " " + "(" + number + ")" + 
                    "\t")
                
                for k, v in process.iteritems():
                    category = k
                    number = v
                    out_file.write(category + " " + "(" + number + ")" + "\t")
                
                for k, v in compartment.iteritems():
                    category = k
                    number = v
                    out_file.write(category + " " + "(" + number + ")" + "\n")

                treatment = line.split("\t")[0]
                description.discard(line.split("\t")[1])
                
                function = {}
                process = {}
                compartment = {}
            
            # Treats the first line of the in_file document
            if count == 1:
                
                treatment = line.split("\t")[0]
                description.discard(line.split("\t")[1])

            contig = line.split("\t")[0]
            blast = line.split("\t")[1]
            
            # Fills in the dictionaries corresponding to the 3 GO categories with the info
            # written in the in_file
            if contig == treatment:

                GO = ""
                term = ""
                
                if line.split("\t")[4] == "F":
                    term = line.split("\t")[3]
                    GO = line.split("\t")[2]
                    function[term] = GO
                    
                if line.split("\t")[4] == "C":
                    term = line.split("\t")[3]
                    GO = line.split("\t")[2]
                    compartment[term] = GO

                if line.split("\t")[4] == "P":
                    term = line.split("\t")[3]
                    GO = line.split("\t")[2]
                    process[term] = GO
