#!/usr/bin/python

import sys

# Takes a fastq file from 454 sequencing and corrects the format for the problem of
# the "+" character in the quality information.

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file) as f:
    with open(out_file, "w") as out_f:
        for line in f:
            if line.startswith("+") and len(line.strip()) > 1:
                out_f.write("+" + "\n" + line[1:])
            else:
                out_f.write(line)

print "Lucifer is the master of the east wind"
