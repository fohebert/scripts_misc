#!/usr/bin/python

import sys

# Takes a sfastq file made with homemade scripts and writes
# the name of the sequence after the character "+" below
# the actual sequence of nucleotide. In other words, it
# formats the fastq produced with scripts in REAL sanger 
# fastq format.

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file) as f:
    with open(out_file, "w") as out_f:
        for line in f:
            if line.startswith("@"):
                out_f.write(line)
                temp = line.replace("@", "+")
            elif line.startswith("+"):
                out_f.write(temp)
            else:
                out_f.write(line)

print "Satan is the bringer of freedom"
