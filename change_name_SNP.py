#!/usr/bin/python

"""

SYNOPSIS: 
    
    Takes the output file from snp_count.py and modifies
    the name of the SNPs which is "contig<number>" for 
    "cocl_<number>_<SNP position>".

USAGE: 
    
    %program <in_file SNP count> <out_file>

CREDITS:

    Doctor Pants 2011 \m/

"""

import sys
import re

try:
    in_file = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(0)

with open(in_file) as f:
    with open(out_file, "w") as f2:
        for line in f:
            
            if line.startswith("Contig_nb"):
                f2.write(line)
            
            else:
            
                # stores the new contig name in the object "SNP_name"
                contig = line.split("\t")[0]
                position = line.split("\t")[1]
                SNP_name = "cocl_" + contig[6:11] + "_" + position
                
                # stores the rest of the information of the line in the object "rest"
                rest = re.split("contig[0-9]*\t", line)[1]
            
                # writes down the new line with the new contig name in the output
                f2.write(SNP_name + "\t" + rest)
