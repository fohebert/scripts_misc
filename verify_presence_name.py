#!/usr/bin/python

"""
SYNOPSIS :
    Takes a text file containing contig or probe
    names and verifies if the name is already
    written in another text file containg many names.
    It writes down the names that are not found in
    the reference file in an output text file.

NB: one name per line for both input files

USAGE : 
    %program <names to verify> <ref names> <output>

Doc Pants 2011 \m/
"""

import sys

try:
    names_to_verify = sys.argv[1] # Text file containing the names that are possibly in another file
    ref_names = sys.argv[2]       # Text file containing reference names - This is the reference file
    output_file = sys.argv[3]     # Output file containing the names that are not in the reference
except:
    print __doc__
    sys.exit(0)

reference_names = set()

with open(ref_names, "r") as f:
    for line in f:
        line = line.strip()
        reference_names.add(line)

with open(names_to_verify, "r") as f:
    with open(output_file, "w") as f2:
        for line in f:
            line = line.strip()
            if line not in reference_names:
                f2.write(line + "\n")
