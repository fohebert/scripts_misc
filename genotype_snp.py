#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

                        User Commands

SYNOPSIS:

    Extracts the genotype of individuals from the output file of
    "eliminate_full_hetero_SNP.py". Generates a matrix of  geno-
    types in a tabulated text file of this type : 
    
    Individual  locus1      locus2      locus3
    indiv_1     A   T       C   G       G   G
    indiv_2     A   A       C   G       G   G
    ...
    
    Input file formats:
    
    SNP_TABLE:
    
    Contig_nb	Pos	tag_name	A	C	G	T	N	*	-
    contig45	11	C-D16	3	0	0	0	0	0	0
    contig45	11	C-D17	1	0	0	2	0	0	0
    contig45	11	C-D18	5	0	0	4	0	0	0
    contig45	11	C-N6	4	0	0	1	0	0	0
    contig45	11	C-N7	4	0	0	1	0	0	0
    contig45	11	C-N8	6	0	0	2	0	0	0
    contig45	11	XX_noTag	1	0	0	0	0	0	0
    
    TAGS:
    
    C-D16
    C-D17
    C-D18
    ...
    
    SNP_NAMES:
    
    The name of each locus (each SNP), one name per line:
    
    cocl_001_11
    cocl_002_35
    cocl_003_44
    ...
    
    NB : use zeros in the numbers in order to be able to
    sort properly the loci and to produce an adequate
    out file.

USAGE:
    
    %program <SNP_table> <tags> <snp_names> <output_file>
    
CREDITS:

    Doctor Pants 2011 \m/
    
"""

# Extracts the genotype of individuals based on the output file from "eliminate_full_hetero_SNP.py"

import sys
import re

try:
    input_snp = sys.argv[1]
    tags = sys.argv[2]
    names_snp = sys.argv[3]
    output = sys.argv[4]
except:
    print __doc__
    sys.exit(0)

# Creates a dictionary of dictionaries : one dictionary per individual. Each dictionary will contain
# the genotype of all the SNPs (all the loci).

dictionaries = {}

with open(tags) as tags:
    for line in tags:
        line = line.strip()
        if line != "":
            dictionaries[line] = {}

# Stores the name of each locus (each SNP) in a list.

snps = []

with open(names_snp) as names:
    with open(output, "w") as out_file:
        
        for line in names:
            line = line.strip()
            if line != "":
                snps.append(line + "\t")

# Writes the genotype of each individual at each locus in a separate dictionary (one dictionary
# per individual).

with open(input_snp) as in_snp:            
    for line in in_snp:        
        
        to_add = ""
        alleles = []
        count = 0
        
        if line.startswith("Contig_nb"):
            continue
        
        if re.search('XX_noTag', line) != None:
            continue
        
        if re.findall('C-[A-Z][0-9]*', line.split("\t")[2]) != "":    
            
            sample = line.split("\t")[2]
            locus = line.split("\t")[0]
                     
            if float(line.split("\t")[3]) > 0:
                alleles.append("A")
                count += 1
            if float(line.split("\t")[4]) > 0:
                alleles.append("C")
                count += 1
            if float(line.split("\t")[5]) > 0:
                alleles.append("G")
                count += 1
            if float(line.split("\t")[6]) > 0:
                alleles.append("T")
                count += 1

            if count == 1:
                alleles.append(alleles[0])
            if count == 0:
                alleles = ["-","-"]
                    
            dictionaries[sample][locus] = alleles
            
# Builds the output file

with open(output, "w") as out_file:
    
    # Writes the word "Individual" in the top left corner of the matrix.
    out_file.write("Individual" + "\t")
    
    # Writes the names of the SNP in the first row.
    for objects in snps:
        out_file.write(objects)
    
    # Write down the genotype information in the out file.
    for sample in list(sorted(dictionaries)):
        
        out_file.write("\n" + sample + "\t")
        
        for locus in list(sorted(dictionaries[sample])):
            out_file.write("/".join(list(sorted(dictionaries[sample][locus]))) + "\t")
