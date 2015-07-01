#!/usr/bin/python

"""

                   User Commands

SYNOPSIS:
    
    Counts the number of SNPs per contig and generates
    an output file containing the contigs names (one
    name per line) and the number of SNPs for this
    contig (tabulated text file).

USAGE:
    
    %program <SNP_file> <contig_names> <output>
    
    SNP file: 
    SNP output from the script "snp_count.py"
    
    contig_names:
    one name per line
    
    output:
    
    contig001   10
    contig002   7
    contig003   13
    ...

CREDITS:
    
    Doctor Pants 2011 \m/

"""

import sys

try:
    SNP_in_file = sys.argv[1]
    contig_names = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

# Creates a set of contig names
names = set()

with open(contig_names) as f:
    for line in f:
        line = line.strip()
        if line != "":
            names.add(line)

# Parse the in_file line by line and counts the SNPs for each
# contig present in the SNP table file (in_file).
treatment = ""
contig = ""
count = 0

with open(SNP_in_file) as in_file:
    with open(out_file, "w") as o_file:
        o_file.write("Contig" + "\t" + "num. SNPs")
        
        for line in in_file:
            
            line = line.strip()
                        
            if line.startswith("Contig_nb"):
                continue
            
            if line.split("\t")[0] in names and count == 0:
                treatment = line.split("\t")[0]
                names.discard(treatment)

            if line.split("\t")[2] == "XX_noTag":
                count += 1
            
            if line.split("\t")[0] != treatment:
                
                o_file.write(treatment + "\t" + str(count) + "\n")
                treatment = line.split("\t")[0]
                names.discard(treatment)
                count = 0
            
            if line.split("\t")[2] == "end":
                o_file.write(treatment + "\t" + str(count) + "\n")
