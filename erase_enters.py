#!/usr/bin/python

import sys
import re

in_file = sys.argv[1]
out_file = sys.argv[2]

with open(in_file, "r") as f:
	with open(out_file, "w") as f2:
		for line in f:
			line = line.strip()
			if line.startswith(">") == True:
				f2.write(line + "\n")
			else:
				f2.write(line)
