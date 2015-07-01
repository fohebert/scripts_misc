#!/usr/bin/python

"""
Takes an ace file and extracts only the wanted contigs\ 
based on a text file containing the name of the sequences.

Usage : %program <ace file> <file with names> <output>
"""

import sys

try:
    ace_file = sys.argv[1]
    wanted_contigs = sys.argv[2]
    output_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

wanted = set()

with open(wanted_contigs, "r") as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

wanted_w = False

with open(ace_file, "r") as f:
    with open(output_file, "w") as f2:
        for line in f:
            line = line.strip()
            if line[0:2] == "AS":
                f2.write(line + "\n")
            if line[0:2] == "CO":
                wanted_w = False
                if line.split()[1] in wanted:
                    wanted.remove(line.split()[1])
                    wanted_w = True
            if wanted_w == True:
                f2.write(line + "\n")
