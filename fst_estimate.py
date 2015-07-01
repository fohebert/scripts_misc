#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                            User Commands

\033[1mDESCRIPTION\033[0m
        This program takes a genotype file created with the script
        "snps-geno_to_bayescan_assembled_genes.py"  and  estimates
        the Fst values for  each locus passing  the threshold of a
        percentage of heterozygous <= 0.65.

\033[1mUSAGE\033[0m
        %program <genotype_file> <num_samples> <output>

\033[1mOUTPUT FILE - FORMAT\033[0m
        One locus by line, tab delimited : locus number + Fst value
        To each locus is  assigned a number that corresponds to the
        same Bayescan input number.

\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys

try:
    in_genotype = sys.argv[1]
    num_samp = int(sys.argv[2])
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(2)

smp_col = {}
loci_count = 1
with open(in_genotype, "rU") as geno_in:
    with open(out_file, "w") as out_f:
        
        out_f.write("SNP num." + "\t" + "Gene" + "\t" + "Gene pos." + "\t" + "Fst" + "\t" + "Hetero tot.")
        
        for line in geno_in:
            line = line.strip()
            
            # Associates the name of the sample to its corresponding column position in a dictionnary
            if line.startswith("Gene"):
                samples = num_samp
                ind_column = 5
                
                while samples > 0:
                    sp = line.split("\t")[ind_column]
                    smp_col[ind_column] = sp
                    ind_column += 1
                    samples -= 1
            
            elif line.startswith("Gene") == False:
                
                a1 = line.split("\t")[2].split("/")[0]
                a2 = line.split("\t")[2].split("/")[-1]
                
                # Loop that calculates for each locus the number of successfully sequenced individuals
                # in total and the number of heterozygous samples amnong this total.
                count = num_samp
                col_num = 5
                a1_tot_counts = 0
                a2_tot_counts = 0
                a1_d_counts = 0
                a2_d_counts = 0
                a1_n_counts = 0
                a2_n_counts = 0
                num_hetero = 0
                tot_sequenced = 0
                while count > 0:
                    if line.split("\t")[col_num] != "0":
                        
                        tot_sequenced += 1
                        
                        # Allele counts for DWARFS
                        if smp_col[col_num].startswith("C-D"):
                            
                            samp_a1 = line.split("\t")[col_num].split("/")[0]
                            samp_a2 = line.split("\t")[col_num].split("/")[1]
                            
                            if samp_a1 != samp_a2:
                                num_hetero += 1
                            
                            # Makes the allele count for the first allele according to the allele identity
                            if samp_a1 == a1:
                                a1_tot_counts += 1
                                a1_d_counts += 1
                            elif samp_a1 == a2:
                                a2_tot_counts += 1
                                a2_d_counts += 1
                            
                            # Makes the allele count for the second allele according to the allele identity
                            if samp_a2 == a1:
                                a1_tot_counts += 1
                                a1_d_counts += 1
                            elif samp_a2 == a2:
                                a2_tot_counts += 1
                                a2_d_counts += 1
                        
                        # Allele counts for NORMALS
                        if smp_col[col_num].startswith("C-N"):
                            
                            samp_a1 = line.split("\t")[col_num].split("/")[0]
                            samp_a2 = line.split("\t")[col_num].split("/")[1]
                            
                            if samp_a1 != samp_a2:
                                num_hetero += 1
                            
                            # Makes the allele count for the first allele according to the allele identity
                            if samp_a1 == a1:
                                a1_tot_counts += 1
                                a1_n_counts += 1
                            elif samp_a1 == a2:
                                a2_tot_counts += 1
                                a2_n_counts += 1
                            
                            # Makes the allele count for the second allele according to the allele identity
                            if samp_a2 == a1:
                                a1_tot_counts += 1
                                a1_n_counts += 1
                            elif samp_a2 == a2:
                                a2_tot_counts += 1
                                a2_n_counts += 1
                    
                    count -= 1
                    col_num += 1
                
                percentage_hetero = float(num_hetero) / float(tot_sequenced)
                
                if percentage_hetero <= 0.65:
                    
                    if loci_count == 23:
                        print smp_col
                        print "a1 = ", a1
                        print "a2 = ", a2
                        print "a1_dwarf = ",a1_d_counts
                        print "a2_dwarf = ",a2_d_counts
                        print "a1_norm = ",a1_n_counts
                        print "a2_norm = ",a2_n_counts
                        print "a1_tot = ", a1_tot_counts
                        print "a2_tot = ", a2_tot_counts
                    
                    # Calculates the allelic frequencies for the dwarfs
                    freq_a1_d = float(a1_d_counts)/float(a1_d_counts + a2_d_counts)
                    freq_a2_d = float(a2_d_counts)/float(a1_d_counts + a2_d_counts)
                    
                    # Calculates the allelics frequencies for the normals
                    freq_a1_n = float(a1_n_counts)/float(a1_n_counts + a2_n_counts)
                    freq_a2_n = float(a2_n_counts)/float(a1_n_counts + a2_n_counts)
                    
                    # Calculates the total allele freqencies
                    freq_a1_tot = float(a1_tot_counts)/float(a1_tot_counts + a2_tot_counts)
                    freq_a2_tot = float(a2_tot_counts)/float(a1_tot_counts + a2_tot_counts)
                    
                    # Calculates the expected heterozygosity in the total population (ht)
                    ht = 2*freq_a1_tot*freq_a2_tot
                    
                    # Calculates the mean expected heterozygosity for the 2 populations (hs)
                    hd = 2*freq_a1_d*freq_a2_d
                    hn = 2*freq_a1_n*freq_a2_n
                    hs = (hd + hn) / float(2)
                    
                    # Calculates the Fst for this locus
                    fst = abs(ht - hs) / ht
                    out_f.write("\n" + str(loci_count) + "\t" + line.split("\t")[0] + "\t" + line.split("\t")[1] + "\t" + str(fst) + "\t" + str(ht))
                
                    loci_count += 1

print "\nDone\n"
