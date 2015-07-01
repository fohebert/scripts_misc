#!/usr/bin/python

# Takes a blast result file parsed with a 'Blast2ID' type script
# and writes in a new text file the names of the queries that
# blasted successfully (leaving behind the "No hits" results)

import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file, "r") as f:
	with open(out_file, "w") as f2:
		for line in f:
			line = line.strip()
			if line.find("No hits found") == -1:
				f2.write(line + "\n")
