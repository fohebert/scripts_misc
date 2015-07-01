#!/usr/bin/python

"""\nExtracts the names of wanted reads from an ACE file

Format of wanted names file : one name per line

Usage : 
	%program <input_ace_file> <wanted_names_file> <name_output_file> 
	
	Doctor Pants 2011 \m/ \n """

import sys
import re

try:
	in_file = sys.argv[1]
	wanted_names = sys.argv[2]
	out_file = sys.argv[3]
except:
	print __doc__
	sys.exit(0)

names = set()

with open(wanted_names) as f:
	for line in f:
		line = line.strip()
		if line != "":
			names.add(line)

write_name = False

with open(in_file) as in_f:
	with open(out_file, "w") as out_f:
		for line in in_f:
			l = line.strip()
			if l.startswith('CO ') and l.split()[1] in names:
				write_name = True
			elif l.startswith('AF C') and write_name == True:
				out_f.write(l.split()[1] + "\n")
			elif l.startswith('CO ') and l.split()[1] not in names:
				write_name = False
