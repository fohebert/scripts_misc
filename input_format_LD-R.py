#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                                       Help Manual

\033[1mDESCRIPTION\033[0m
        Convert a genotype matrix file in the  appropriate input  format for R in order 
        to perform a LD analysis.

\033[1mUSAGE\033[0m
        %program <genotype_matrix> <perc_hetero> <num_samples> <out_dwarf> <out_normal>
            <out_loci_list>
"""

import sys

try:
    genotype_matrix = sys.argv[1]
    perc_hetero = float(sys.argv[2])
    num_smp = int(sys.argv[3])
    output_dwarf = sys.argv[4]
    output_normal = sys.argv[5]
    out_loci = sys.argv[6]
except:
    print __doc__
    sys.exit(2)

smp_col = {}
locus_num = 1
loci = {}
list_loci = []
with open(genotype_matrix, "rU") as geno_mat:
    with open(output_dwarf, "w") as out_d:
        with open(output_normal, "w") as out_n:
            for line in geno_mat:
                line = line.strip()
                
                if line.startswith("Gene"):
                    
                    num_ind = num_smp
                    ind_col = 5
                    while num_ind > 0:
                        sp = line.split("\t")[ind_col]
                        smp_col[ind_col] = sp
                        ind_col += 1
                        num_ind -= 1
                    
                if line.startswith("Gene") == False:
                    
                    gene = line.split("\t")[0]
                    pos = int(line.split("\t")[1])
                    
                    num_ind = num_smp
                    ind_col = 5
                    num_h_d = 0
                    num_d = 0
                    num_h_n = 0
                    num_n = 0
                    while num_ind > 0:
                        
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
                        num_ind -= 1
                    
                    tot_hetero = float(num_h_d + num_h_n)/float(num_d + num_n)
                    if tot_hetero <= perc_hetero:
                        
                        locus = "loc" + str("%04i" % locus_num)
                        list_loci.append(locus)
                        locus_num += 1
                        
                        if gene not in list(sorted(loci)):
                            
                            loci[gene] = {}
                            loci[gene][pos] = locus
                        
                        elif gene in list(sorted(loci)):
                            
                            loci[gene][pos] = locus
                        
                        num_ind = num_smp
                        ind_col = 5
                        geno_dwf = []
                        geno_nml = []
                        while num_ind > 0:
                            
                            smp = smp_col[ind_col]
                            if line.split("\t")[ind_col] != "0":
                                
                                if smp.startswith("C-D"):
                                    geno = "'" + line.split("\t")[ind_col] + "'"
                                    geno_dwf.append(geno)
                                
                                if smp.startswith("C-N"):
                                    geno = "'" + line.split("\t")[ind_col] + "'"
                                    geno_nml.append(geno)
                                    
                            if line.split("\t")[ind_col] == "0":
                                
                                geno = "NA"
                                if smp.startswith("C-D"):
                                    geno_dwf.append(geno)
                                
                                if smp.startswith("C-N"):
                                    geno_nml.append(geno)
                            
                            ind_col += 1
                            num_ind -= 1
                        
                        out_d.write(locus + " <- genotype(c(" + ",".join(geno_dwf) + "))" + "\n")
                        out_n.write(locus + " <- genotype(c(" + ",".join(geno_nml) + "))" + "\n")
                        
            out_d.write("data <- makeGenotypes(data.frame(" + ",".join(sorted(list_loci)) + "))")
            out_n.write("data <- makeGenotypes(data.frame(" + ",".join(sorted(list_loci)) + "))")
            
            with open(out_loci, "w") as out_l:
                for gene in list(sorted(loci)):
                    for pos in sorted(loci[gene]):
                        out_l.write(loci[gene][pos] + "\t" + gene + "\t" + str(pos) + "\n")
                        
print "\n\033[1mDone!\033[0m\n"
