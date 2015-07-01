#!/usr/bin/python

"""

DESCRIPTION
        Takes a tab text file containing information on 
        the coverage per position (format : position + 
        tab + coverage) and formats it to have an input
        file for R in order to produce a histogram plotting
        the coverage per position.

USAGE
        %program <raw_format> <result_file>
        
        raw_format
        ****************
        Format as outputed by the command line using samtools
        that is used to generate a text file containing the
        positions and their coverage (one position per line)
        from the .bam format produced by CLC.
        
        ...
        145  12
        146  20
        147  19
        148  18
        149  21
        ...
        
        result_file
        ****************
        One column, each line containing a number corresponding
        to the position on the contig and this position is 
        repeated X number of times, X = coverage.

CREDITS
        Doctor Pants 2011 \m/

"""

import sys

try:
    raw_format = sys.argv[1]
    histogram_formated = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

with open(raw_format) as raw:
    with open(histogram_formated, "w") as formated:
        for line in raw:
            
            count = 0
            line = line.strip()
            
            position = line.split("\t")[0]
            coverage = int(line.split("\t")[1])
            
            while coverage > 0:
                coverage -= 1
                formated.write(position + "\n")
