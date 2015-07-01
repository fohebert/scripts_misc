#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

                        User Commands

SYNOPSIS:

    Takes the output file of the script "snp_count.py" and
    eliminates the loci that only have heterozygous individuals.

USAGE:
    
    %program <SNP_table from "snp_count.py"> <output>

CREDITS:

    Doctor Pants 2011 \m/
    
"""

import sys
import re

try:
    input_snp = sys.argv[1]
    output = sys.argv[2]

except:
    print __doc__
    sys.exit(0)

end = False
discard = False
heterozygous = 0
allele_count = 0
SNP = []

with open(input_snp, "r") as in_file:
    with open(output, "w") as out_file:
        for line in in_file:

            if line.startswith("Contig_nb"):
                out_file.write(line)

            else:
                            
                SNP.append(line)
                
                allele_count = 0
                end = False
                discard = False
                                            
                if re.findall('[A-Z]-[A-Z][0-9]*', line.split("\t")[2]) != "":

                    if float(line.split("\t")[3]) != 0:
                        allele_count += 1
                    if float(line.split("\t")[4]) != 0:
                        allele_count += 1
                    if float(line.split("\t")[5]) != 0:
                        allele_count += 1
                    if float(line.split("\t")[6]) != 0:
                        allele_count += 1
                    if allele_count > 1:
                        heterozygous += 1
            
                if line.split("\t")[2] == 'XX_noTag':
                    end = True
                    
                    if heterozygous == 6:
                        discard = True

            if end:
                if discard:
                    SNP = []
                else:
                    out_file.write(SNP[0]+SNP[1]+SNP[2]+SNP[3]+SNP[4]+SNP[5]+SNP[6])
                    SNP = []
                heterozygous = 0
