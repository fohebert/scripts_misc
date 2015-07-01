#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
                            User Commands

\033[1mDESCRIPTION\033[0m
        Takes the genotype matrix file (tab delimited) generated with
        the script "snps-geno_to_bayescan_assembled_genes.py" and re-
        turns  in an output file the heterozygosity for each SNP that
        has Ho <= 0.65 PER POPULATION.

\033[1mUSAGE\033[0m
        %program <genotypes> <num_samples> <output>

\033[1mOUTPUT FILE - DESCRIPTION\033[0m
        
        Gene      pos      Ho dwarf     Ho Normal
        HSP-90    435      0.1662       0.562
        ...
        
\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys

try:
    genotype_file = sys.argv[1]
    num_samples = int(sys.argv[2])
    output_file = sys.argv[3]
except:
    print __doc__
    sys.exit(2)

smp_col = {}
hetero_d = 0
hetero_n = 0
with open(genotype_file, "rU") as geno_f:
    with open(output_file, "w") as out_f:
        
        out_f.write("Gene" + "\t" + "pos" + "\t" + "Ho dwarf" + "\t" + "Ho normal")
        
        for line in geno_f:
            line = line.strip()
            
            if line.startswith("Gene"):
                samples = num_samples
                ind_column = 5
                
                while samples > 0:
                    sp = line.split("\t")[ind_column]
                    smp_col[ind_column] = sp
                    ind_column += 1
                    samples -= 1
            
            if line.startswith("Gene") == False:
                
                gene = line.split("\t")[0]
                pos = line.split("\t")[1]
                
                ind_col = 5
                num_smp = num_samples
                num_h_d = 0
                num_d = 0
                num_h_n = 0
                num_n = 0
                while num_smp > 0:
                    
                    if line.split("\t")[ind_col] != "0":
                        
                        smp = smp_col[ind_col]
                        if smp.startswith("C-D"):
                            a1 = line.split("\t")[ind_col].split("/")[0]
                            a2 = line.split("\t")[ind_col].split("/")[-1]
                            if a1 != a2:
                                num_h_d += 1
                            num_d += 1
                        
                        if smp.startswith("C-N"):
                            a1 = line.split("\t")[ind_col].split("/")[0]
                            a2 = line.split("\t")[ind_col].split("/")[-1]
                            if a1 != a2:
                                num_h_n += 1
                            num_n += 1
                    
                    ind_col += 1
                    num_smp -= 1
                
                Ho_d = float(num_h_d)/float(num_d)
                Ho_n = float(num_h_n)/float(num_n)
                Ho_tot = float(num_h_d + num_h_n)/float(num_n + num_d)
                
                if Ho_tot <= 0.65:
                    out_f.write("\n" + gene + "\t" + pos + "\t" + str(Ho_d) + "\t" + str(Ho_n))

print "\n\033[1mDone !\033[0m\n"
