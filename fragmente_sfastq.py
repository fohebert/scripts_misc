#!/usr/bin/python

import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

count = 0

with open(in_file, "r") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            line = line.strip()
            if line != "":
                if line.startswith("@C-"):
                    count += 1
                    if count <= 500000:
                        out_f.write(line + "\n")
                else:
                    if count <= 500000:
                        out_f.write(line + "\n")
