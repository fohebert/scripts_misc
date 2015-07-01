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
        %program <blast_output> <e-value threshold (int)> <output>

\033[1mINPUT FILE - FORMAT\033[0m
        
        Blast output
        ************
        It has to be either the blast result file (.txt) from NCBI's  website  OR
        the blast result file from BLAST+ (e.g local  blast,  swissprot,  nr,  nt
        or custom database) in FORMAT 0 (default).
        
        The program will look for the best blast result to output the exon-intron
        breaks, though it is possible to give it an  input file in which ALL  the
        blast results are good and must be  kept (e.g if a set  of contigs or  se-
        quences from an exon capture chip have been blasted against the target se-
        quences). The e-value to give as input information must follow the format:
        
        [0-9]*e-[0-9]*
        
        e.g : 6e-38

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
    evalue = float(sys.argv[2])
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(2)

contigs = {}
contig = ""
target = False
started = False

with open(blast_in, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        if line.startswith("Query="):
            contig = line.split("Query= ")[1].split("_")[0] + "_" + str("%05i" % \
                int(line.split("Query= ")[1].split("_")[1]))
            contigs[contig] = set()
            started = False
            target = False
        
        elif line.startswith("> "):
            if target:
                started = False
            elif target == False:
                target = True
        
        elif re.findall("Expect =", line) != [] and target:
            current_eval = float(line.split()[-1])
            if current_eval <= evalue:
                started = True
        
        elif re.findall("Query  [0-9]*", line) != [] and started:
            pos1 = int(re.split(" [0-9]* ", line)[1])
            pos2 = int(re.split(" [0-9]* ", line)[-1])
            
            for pos in range(pos1,pos2+1):
                contigs[contig].add(pos)

exon = set()
with open(out_file, "w") as out_f:
    
    out_f.write("Contig" + "\t" + "Exon positions")
    
    for contig in list(sorted(contigs)):
        out_f.write("\n" + contig + "\t")
        
        lgth = len(contigs[contig])
        first_pos = 0
        for i in range(0,lgth):
            
            if i == 0:
                first_pos = sorted(contigs[contig])[i]
            
            elif sorted(contigs[contig])[i] == sorted(contigs[contig])[i-1] + 1:
                pass
            
            elif sorted(contigs[contig])[i] > sorted(contigs[contig])[i-1] + 1:
                out_f.write(str(first_pos) + "-" + str(sorted(contigs[contig])[i-1]) + ";")
                first_pos = sorted(contigs[contig])[i]
            
            if i == (lgth - 1):
                out_f.write(str(first_pos) + "-" + str(sorted(contigs[contig])[i]) + ";")
                
    print "\nFile processed\n"
