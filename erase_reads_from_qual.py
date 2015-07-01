#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

"""\nTakes a qual file and extracts all the sequences except\
    the ones specified  in a text file, one name per line. 
    
    Usage : <qual file> <unwanted names> <output file>\n"""

try:
    qual_file = sys.argv[1]
    unwanted_names = sys.argv[2]
    out_file = sys.argv[3]
    
except:
    print __doc__
    sys.exit(0)

unwanted = set()
with open(unwanted_names) as f:
    for line in f:
        line = line.strip()
        if line != "":
            unwanted.add(line)

qual = SeqIO.parse(open(qual_file), "qual")
end = False
with open(out_file, "w") as f2:
    for seq in qual:
        if seq.id not in unwanted:
            SeqIO.write([seq], f2, "qual")
    print '\nDone !\n'
