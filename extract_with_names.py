#!/usr/bin/python

import sys
import re
from Bio import SeqIO

# Takes a text file containing many blasted contigs (from a "de nevo" assembly) and extracts the name 
# of the desired ones based on their gene ID.

contig_file = sys.argv[1] # Tab delimited text file containing the names of the contigs with gene IDs per contig (from NCBI)
text_wanted = sys.argv[2] # Text file containing a list of the wanted gene ID
output_file = sys.argv[3] # Output file containing the names of the wanted contigs with their gene IDs. One contig per line.

wanted = set()
final = []

with open(text_wanted) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

with open(output_file, "w") as outf:
    with open(contig_file) as f2:
        for line in f2:
            line = line.strip()
            for gene_id in wanted:
                if line.find(gene_id) != -1:
                    outf.write(line + "\n")
