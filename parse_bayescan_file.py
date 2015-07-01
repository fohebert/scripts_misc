#!/usr/bin/python

"""

Parses a bayescan input file and discard the loci for which
not all individuals have been genotyped. Writes down only the
loci that are fully genotyped for all samples in an output file
of the same format as the input file.

Usage: %program <input_bayescan format> <unwanted loci> <output>

unwanted loci = list of the unwanted loci (one loci number per line)

IMPORTANT : CHANGE THE NUMBER OF LOCI IN THE FIRST LINE OF THE
OUPUT FILE !!!!!!!!!!!!!!!

"""

import sys
import re

try:
    in_file = sys.argv[1]
    list_file = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

unwanted_loci = set()
with open(list_file,"r") as in_f:
    for line in in_f:
        line = line.strip()
        if line != "":
            unwanted_loci.add(line)

with open(in_file, "r") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            
            if line[0:6] == "[loci]":
                out_f.write(line + "\n")
            if line[0:6] == "[Popul":
                out_f.write(line)
            if line[0:6] == "[pop]=":
                out_f.write("\n" + line + "\n")
            
            if re.findall("^[1-9]", line):
                if line.split()[0] not in unwanted_loci:
                    out_f.write(line)
