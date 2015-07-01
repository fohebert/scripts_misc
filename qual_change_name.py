#!/usr/bin/env python

"""
\033[1mDESCRIPTION\033[0m
    Takes a .QUAL file and changes the name when there is
    multiple words with spaces. Keeps only the first word
    and  discards  everything  else that  is separated by
    either a tab or a space.

\033[1mUSAGE\033[0m
    %program <input_qual> <output_qual>
"""

import sys

try:
    in_qual = sys.argv[1]
    out_qual = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

count = 0
with open(in_qual, "rU") as in_q:
    with open(out_qual, "w") as out_q:
        for line in in_q:
            line = line.strip()
            
            if line.startswith(">"):
                name = line.split(">")[-1].split()[0]
                if count == 0:
                    out_q.write(">" + name)
                else:
                    out_q.write("\n" + ">" + name)
                count +=1 
            else:
                out_q.write("\n" + line)

print "\n\033[1mJob Done\033[0m\n"
