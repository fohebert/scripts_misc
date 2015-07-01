#!/usr/bin/python

"""
\033[1mDESCRIPTION\033[0m
    Takes a text file with one name (or whatever info on the line)
    per line and transposes into columns. In other words, it
    transposes rows into columns.

\033[1mUSAGE\033[0m
    %program <in_rows> <out_columns>
"""

import sys

try:
    in_file = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

lst_sp = []
with open(in_file, "rU") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            line = line.strip()
            lst_sp.append(line)
        out_f.write("\t" + "\t".join(lst_sp) + "\n")
