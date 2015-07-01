#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

                        User Commands

\033[1mUSAGE\033[0m:
    
    %program <SNP_count> <tags> <loci> <CLC_SNP_table> <output_file>

\033[1mSYNOPSIS\033[0m:

    Extracts the genotype of individuals from the 
    output of the script "SNP_count.py" and converts
    it into genepop format. Then, the genepop format
    can be imported in the programs BAYESFST or PGD
    Spider to get converted into another format (e.g.
    BAYESCAN).
    
    Input file formats:
    
    \033[1mSNP_count\033[0m:
    
    Contig_nb	Pos	tag_name	A	C	G	T	N	*	-
    contig45	11	C-D16	3	0	0	0	0	0	0
    contig45	11	C-D17	1	0	0	2	0	0	0
    contig45	11	C-D18	5	0	0	4	0	0	0
    contig45	11	C-N6	4	0	0	1	0	0	0
    contig45	11	C-N7	4	0	0	1	0	0	0
    contig45	11	C-N8	6	0	0	2	0	0	0
    contig45	11	XX_noTag	1	0	0	0	0	0	0
    
    \033[1mTAGS\033[0m:
    
    C-D16
    C-D17
    C-D18
    ...
    
    \033[1mLOCI\033[0m:
    
    The name of each locus (each SNP), one name per line:
    
    cocl_001_11
    cocl_002_35
    cocl_003_44
    ...
    
    NB : use zeros in the numbers in order to be able to
    sort properly the loci and to produce an adequate
    out file.
    
    \033[1mCLC_SNP_table\033[0m:
    
    CLC SNP table parsed for the desired SNPs, i.e only
    the ones that are potentially interesting. It is pos-
    sible to use the script "parse_CLC_SNP_table.py" to
    obtain this file.
    
    MAKE SURE THE NAMES ARE ALL THE SAME :
    
       --> cocl_contig#_position <--
    
    \033[1mOutput file format\033[0m:
    
    ________the file starts below this line___________
    Title line : 'Exon capture chip, SNPs - Whitefish'
    cocl_001_11 , cocl_002_35 , cocl_003_44
    POP
    C-D16  0101 0102 0000
    C-D17  0102 0102 0102
    C-D18  0101 0101 0202
    POP
    C-N6   0102 0202 0101
    C-N7   0202 0202 0202
    C-N8   0102 0102 0102
    _________the file ends above this line____________
    
\033[1mCREDITS\033[0m:

    Doctor Pants 2011 \m/
    
"""

# Extracts the genotype of individuals based on the output file from "eliminate_full_hetero_SNP.py"

import sys
import re

try:
    input_snp = sys.argv[1]
    tags = sys.argv[2]
    names_snp = sys.argv[3]
    CLC_SNP_table = sys.argv[4]
    output = sys.argv[5]
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

# Stores the genotype information at each locus

locus_genotype = {}

with open(CLC_SNP_table, "r") as in_f:
    for line in in_f:
        line = line.strip()
        if line != "":
            locus_n = "cocl_" + line.split("g")[1].split()[0] + "_" + line.split()[1]
            genotype = line.split()[5] + line.split()[14]
            locus_genotype[locus_n] = [genotype[0],genotype[1]]

# Writes the genotype of each individual at each locus in a separate dictionary (one dictionary
# per individual).

with open(input_snp) as in_snp:            
    for line in in_snp:        
        
        to_add = ""
        alleles = []
        genotype = []
        count = 0

        A1 = ""
        A2 = ""
        
        if line.startswith("Contig_nb"):
            continue
        
        if re.search('XX_noTag', line) != None:
            continue
        
        if re.findall('C-[A-Z][0-9]*', line.split("\t")[2]) != "":    
            
            sample = line.split("\t")[2]
            locus = line.split("\t")[0]
                     
            if float(line.split("\t")[3]) > 0:
                genotype.append("A")
                count += 1
            if float(line.split("\t")[4]) > 0:
                genotype.append("C")
                count += 1
            if float(line.split("\t")[5]) > 0:
                genotype.append("G")
                count += 1
            if float(line.split("\t")[6]) > 0:
                genotype.append("T")
                count += 1
                
            if count == 0:
                genotype = ["00","00"]
                alleles = genotype
                dictionaries[sample][locus] = alleles
            else:
                if count == 1:
                    genotype.append(genotype[0])
                
                if genotype[0] == locus_genotype[locus][0]:
                    A1 = "01"
                if genotype[0] == locus_genotype[locus][1]:
                    A1 = "02"
                if genotype[1] == locus_genotype[locus][0]:
                    A2 = "01"
                if genotype[1] == locus_genotype[locus][1]:
                    A2 = "02"
            
                alleles.append(A1)
                alleles.append(A2)
                
                dictionaries[sample][locus] = alleles
            
# Builds the output file

dwarf = 0
normal = 0

with open(output, "w") as out_file:
    
    # Writes the word "Individual" in the top left corner of the matrix.
    out_file.write("Title line : 'Exon capture chip, SNPs - Whitefish'" + "\n")
    
    # Writes the loci on the second line of the output file (separated by comas).
    for objects in snps:
        out_file.write(" " + objects + ",")
    
    # Write down the genotype information in the out file.
    out_file.write("\n" + "POP" + "\n")
    
    for sample in list(sorted(dictionaries)):

        # Writes down the DWARF genotypes
        if sample[2] == "D":
            
            dwarf += 1
            
            if dwarf == 1:
                out_file.write(sample + " " + "," + " ")
            else:
                out_file.write("\n" + sample + " " + "," + " ")
            
            for locus in list(sorted(dictionaries[sample])):
                
                allele1 = dictionaries[sample][locus][0]
                allele2 = dictionaries[sample][locus][1]
                
                out_file.write(allele1 + allele2 + " ")
        
        # Writes down the NORMAL genotypes
        if sample[2] == "N":
            
            normal += 1
            
            if normal == 1:
                out_file.write("\n" + "POP" + "\n" + sample + " " + "," + " ")
            else:
                out_file.write("\n" + sample + " " + "," + " ")
                
            for locus in list(sorted(dictionaries[sample])):
                allele1 = dictionaries[sample][locus][0]
                allele2 = dictionaries[sample][locus][1]
                
                out_file.write(allele1 + allele2 + " ")
    
    print "\nSatan represents vital existence instead of spiritual pipe dreams\n"
