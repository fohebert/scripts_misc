#!/usr/bin/python

"""
                        User Commands

\033[1mSYNOPSIS\033[0m:
        Takes a file and looks line per line to see if there
        is any lowercases. If there is lowercases, it changes
        them for uppercases and creates the exact same file
        as input file in an output file, but all the letters
        are uppercase instead of lowercase.
        
\033[1mUSAGE\033[0m:
        $program <input_lowercase_file> <output_uppercase_file>
        
\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys

try:
    lowercase_file = sys.argv[1]
    output_uppercase = sys.argv[2]
except:
    print __doc__
    sys.exit(1)

with open(lowercase_file, "r") as in_f:
    with open(output_uppercase, "w") as out_f:
        for line in in_f:
            line = line.strip()
            if line != "":
                line = line.upper()
                out_f.write(line + "\n")
