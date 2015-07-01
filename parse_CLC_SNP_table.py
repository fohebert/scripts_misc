#!/usr/bin/python

"""

                           User Commands

\033[1mDESCRIPTION\033[0m
        Parses the output file "SNP_table" from CLC and generates
        a file containing only the wanted SNPs. The output file 
        is exactly similar to the input file, but it contains on-
        ly the significant SNPs.
        
\033[1mUSAGE\033[0m
        %program <names> <reference> <result_file>
        
        <names> (contig + position, tab text)
        *************************************
        contig00045  11
        contig00045  19
        contig00045  35
        contig00123  141
        ...
        
        <reference>
        *************************************
        contig00045  11  11  SNP  1  A  2  A/T  71,0/29,0  22/9  ...
        contig00045  19  19  SNP  1  A  2  A/T  75,0/25,0  24/8  ...
        contig00045  21  21  SNP  1  G  2  G/A  73,0/27,0  27/10 ...
        contig00045  41  41  SNP  1  C  2  C/T  71,1/28,9  27/11 ...
        ...
        
        <output> 
        *************************************
        identical to <reference>
        
\033[1mCREDITS\033[0m
        Doctor Pants 2011 \m/
        
"""

import sys

try:
    names = sys.argv[1]
    reference = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

contig_names = set()

with open(names) as f:
    for line in f:
        line = line.strip()
        if line != "":
            contig_names.add(line)

with open(reference) as f:
    with open(out_file, "w") as f2:
        line = line.strip()
        for line in f:
            contig_pos = line.split("\t")[0] + "\t" + line.split("\t")[1]
            if contig_pos in contig_names:
                f2.write(line)
