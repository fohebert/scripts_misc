#!/usr/bin/python

import sys
import re

# Takes a text file containing many blasted contigs (from a "de nevo" assembly) and extracts the undesired
# ones (e.g contigs associated with retransposons, rDNA, mtDNA, etc.) based on a text file containing the names
# of the unwanted genes.

contig_file = sys.argv[1] # Tab delimited text file containing the names of the contigs with gene IDs per contig (from NCBI)
text_unwanted = sys.argv[2] # Text file containing a list of the unwanted genes
output_file = sys.argv[3] # Output file containing the names of the unwanted contigs with their gene IDs. One contig per line.

unwanted = set()
final = []

with open(text_unwanted) as f:
    for line in f:
        line = line.strip()
        if line != "":
            unwanted.add(line)

with open(output_file, "w") as outf:
    with open(contig_file) as f2:
        for line in f2:
            line = line.strip()
            for gene_id in unwanted:
                if line.find(gene_id) != -1:
                    outf.write(line + "\n")
