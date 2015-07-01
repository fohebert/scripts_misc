#!/usr/bin/python

# Takes a text file containing 3 tab separated columns : contig name + length + 
# average coverage. One contig per line. Keeps the ones that have a coverage
# superior to 10.

import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file, "r") as f:
	with open(out_file, "w") as f2:
		for line in f:
			l = line.strip()
			name, length, cov = l.split()
			if float(cov) > 10:
				f2.write(name + "\t" + cov + "\n")
