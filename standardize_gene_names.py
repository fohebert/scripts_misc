#!/usr/bin/env python

"""
                            User Commands

SYNOPSIS
        This program parses a gene IDs file to standardise the names
        of the genes corresponding to each contig analysed. The input
        file is one contig per line + tab + gene ID + gene name.
        
        The program eliminates the useless information in the gene names
        (information between brackets relating to species name or other
        useless stuff).

USAGE
        %program <input_gene-IDs> <output_standardised>

CREDITS
        Doctor Pants 2012 \m/
"""

import sys
import re

try:
    gene_IDs = sys.argv[1]
    output_standardised = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

with open(gene_IDs, "rU") as in_f:
    with open(output_standardised, "w") as out_f:
        for line in in_f:
            line = line.strip()
            contig = line.split("\t")[0]
            
            if line.find("; Alt") != -1:
                gene_name = line.split("\t")[1].split(";")[0]
                out_f.write(contig + "\t" + gene_name + "\n")
            
            elif line.find("; Short") != -1:
                gene_name = line.split("\t")[1].split(";")[0]
                out_f.write(contig + "\t" + gene_name + "\n")
            
            elif line.find(", mRNA") != -1:
                gene_name = line.split("\t")[1].split(", mRNA")[0]
                out_f.write(contig + "\t" + gene_name + "\n")
            
            elif re.findall("clone [a-z]*-[a-z]*-[0-9]*-[0-9]*, novel cds", line) != []:
                gene_name = re.split(", ", line)[1]
                out_f.write(contig + "\t" + gene_name + "\n")
            
            elif re.findall("clone [a-z]*-[a-z]*-[0-9]*-[0-9]*", line) != []:
                gene_name = re.split("clone [a-z]*-[a-z]*-[0-9]*-[0-9]* ", line)[1]
                out_f.write(contig + "\t" + gene_name + "\n")
            
            else:
                if line != "":
                    gene_name = line.split("\t")[1]
                    out_f.write(contig + "\t" + gene_name + "\n")
