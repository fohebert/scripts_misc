#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
                                  User Commands

\033[1mSYNOPSIS\033[0m

        Formats the file from the script "VCF_parser.py" to create an input file
        for the program BAYESCAN.

\033[1mUSAGE\033[0m

        %program <SNPs_genotypes> <sample_num>  <max_percent_hetero>  <SNP_info>
        <bayescan>
        
        NB : Percentage of heterozygous allowed to keep for the genotype count. It
        is a FLOAT number between 0.0 and 1.0
        
\033[1mOUTPUT FILE - FORMAT\033[0m
        
        file: out_freq
        ***************
        
        __________file starts below the line__________
        
        Gene      Pos.  num.  Geno.  freq_a1  freq_a2   all_freq_div   p-value
        gene_x    945   1      A/G    0.20    0.92        0.72         0.03264
        gene_y    ...
        
        ___________file ends above the line___________
        
        file: bayescan
        **************

        __________file starts below the line__________
        
        [loci]=508

        [Populations]=2

        [pop]=1
        1  6  2  2 4 
        2  6  2  2 4 
        3  6  2  2 4 
        4  6  2  2 4 
        ...
        
        [pop]=2
        1  6  2  3 3 
        2  6  2  3 3 
        3  6  2  3 3 
        4  6  2  2 4 
        ... 
        
        ___________file ends above the line___________

\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys

try:
    from fisher import pvalue
except:
    print "ERROR - This program needs the package 'Fisher-0.1.4' : http://pypi.python.org/pypi/fisher/"
    sys.exit(2)

try:
    in_file = sys.argv[1]
    num_samp = int(sys.argv[2])
    max_percent_hetero = float(sys.argv[3])
    SNP_info = sys.argv[4]
    out_file = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

# Calculates the pvalue for a 2x2 contingency table using the allelic frequencies in both populations
# Col. = alleles, rows = populations. ***FISHER'S EXACT TEST***
def pval(d1, d2, n1, n2):
    val = 0
    mat = [[d1,d2], [n1,n2]]
    p = pvalue(d1,d2,n1,n2)
    val += p.two_tail
    return val

loci_pop1 = {}
loci_pop2 = {}
smp_col = {}
loci_num = 0
all_percent_hetero = []
with open(in_file, "rU") as in_f:
    with open(SNP_info, "w") as info:
        # Creates the file header for SNP info
        info.write("Gene" + "\t" + "SNP pos." + "\t" + "SNP num." + "\t" + "Geno." + "\t" + "freq_A1_D" + "\t" + 
            "freq_A1_N" + "\t" + "all_freq_div." + "\t" + "p-value")
        
        for line in in_f:
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
                
                # Loop that calculates for each locus the number of successfully sequenced individuals
                # in total and the number of heterozygous samples amnong this total.
                count = num_samp
                col_num = 5
                num_hetero = 0
                num_sequenced = 0
                while count > 0:
                    if line.split("\t")[col_num] != "0":
                        num_sequenced += 1
                        samp_a1 = line.split("\t")[col_num].split("/")[0]
                        samp_a2 = line.split("\t")[col_num].split("/")[1]
                        
                        if samp_a1 != samp_a2:
                            num_hetero += 1

                    count -= 1
                    col_num += 1
                
                # Calculates the percentage of heterozygous individuals throughout both populations
                percent_hetero = float(num_hetero)/float(num_sequenced)
                all_percent_hetero.append(percent_hetero)
                
                if percent_hetero <= max_percent_hetero:
                    
                    gene = line.split("\t")[0]
                    snp_pos = line.split("\t")[1]
                    geno = line.split("\t")[2]
                    
                    # Calculates the allelic frequencies for the REF allele for both populations
                    freq_a1_D = float(line.split("\t")[3].split("/")[0])/float(int(line.split("\t")[3].split("/")[0]) + 
                        int(line.split("\t")[3].split("/")[1]))
                    freq_a1_N = float(line.split("\t")[4].split("/")[0])/float(int(line.split("\t")[4].split("/")[0]) + 
                        int(line.split("\t")[4].split("/")[1]))
                    
                    # Calculates the allelic frequency divergence
                    all_freq_div = abs(freq_a1_D - freq_a1_N)
                    
                    loci_num += 1
                    loci_pop1[loci_num] = {}
                    loci_pop2[loci_num] = {}
                    
                    # Fills in the dictionnary containing the info. on each SNP. Will be used to format the output file
                    # that will be the input file for Bayescan. Every SNP gets attributed a number (from 1 to num_SNPs)
                    # and for each SNP in each population, the allele counts are stored in the dictionnary.
                    tot_num_genes_pop1 = int(line.split("\t")[3].split("/")[0]) + int(line.split("\t")[3].split("/")[1])
                    tot_num_genes_pop2 = int(line.split("\t")[4].split("/")[0]) + int(line.split("\t")[4].split("/")[1])
                    
                    loci_pop1[loci_num]["num_genes"] = 0
                    loci_pop1[loci_num]["num_genes"] += tot_num_genes_pop1
                    loci_pop2[loci_num]["num_genes"] = 0
                    loci_pop2[loci_num]["num_genes"] += tot_num_genes_pop2
                    
                    loci_pop1[loci_num]["a1"] = 0
                    loci_pop1[loci_num]["a1"] += int(line.split("\t")[3].split("/")[0])
                    D_a1 = int(line.split("\t")[3].split("/")[0])
                    loci_pop1[loci_num]["a2"] = 0
                    loci_pop1[loci_num]["a2"] += int(line.split("\t")[3].split("/")[1])
                    D_a2 = int(line.split("\t")[3].split("/")[1])
                    
                    loci_pop2[loci_num]["a1"] = 0
                    loci_pop2[loci_num]["a1"] += int(line.split("\t")[4].split("/")[0])
                    N_a1 = int(line.split("\t")[4].split("/")[0])
                    loci_pop2[loci_num]["a2"] = 0
                    loci_pop2[loci_num]["a2"] += int(line.split("\t")[4].split("/")[1])
                    N_a2 = int(line.split("\t")[4].split("/")[1])
                    
                    # p-value for a Fisher's exact test. Verifies if there is a statistically significant difference in 
                    # the distribution of the allelic frequencies among both populations
                    p_val = pval(D_a1, D_a2, N_a1, N_a2)
                    
                    # Writes the information on the allele frequencies in the output file reserved for that.
                    info.write("\n" + gene + "\t" + snp_pos + "\t" + str(loci_num) + "\t" + geno + "\t" + str(freq_a1_D) + "\t" + 
                        str(freq_a1_N) + "\t" + str(all_freq_div) + "\t" + str(p_val))
        
        global_percent_hetero = sum(all_percent_hetero)/float(len(all_percent_hetero))
        print "\nGlobal percentage of heterozygous among ALL loci: ", global_percent_hetero, "\n"

# Creates the file that will be the input file for Bayescan according to the information stored in the dictionnaries "loci_pop1"
# and "loci_pop2"
with open(out_file, "w") as out_f:
    out_f.write("[loci]=" + str(len(loci_pop1)) + "\n" + "\n" + "[populations]=2" + "\n" + "\n" + "[pop]=1")
    
    for locus in sorted(loci_pop1):
        out_f.write("\n")
        
        if locus < 10:
            out_f.write("   " + str(locus) + "  " + str(loci_pop1[locus]["num_genes"]) + "  " + "2")
        elif 10 <= locus < 100:
            out_f.write("  " + str(locus) + "  " + str(loci_pop1[locus]["num_genes"]) + "  " + "2")
        elif 100 <= locus < 1000:
            out_f.write(" " + str(locus) + "  " + str(loci_pop1[locus]["num_genes"]) + "  " + "2")
        elif locus >= 1000:
            out_f.write(str(locus) + "  " + str(loci_pop1[locus]["num_genes"]) + "  " + "2")
        
        if loci_pop1[locus]["a1"] < 10:
            out_f.write("    " + str(loci_pop1[locus]["a1"]))
        elif loci_pop1[locus]["a1"] >= 10:
            out_f.write("   " + str(loci_pop1[locus]["a1"]))
        
        out_f.write("  " + str(loci_pop1[locus]["a2"]))
        
    out_f.write("\n" + "\n" + "[pop]=2")
    
    for locus in sorted(loci_pop2):
        out_f.write("\n")
        
        if locus < 10:
            out_f.write("   " + str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        elif 10 <= locus < 100:
            out_f.write("  " + str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        elif 100 <= locus < 1000:
            out_f.write(" " + str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        elif locus >= 1000:
            out_f.write(str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        
        if loci_pop2[locus]["a1"] < 10:
            out_f.write("    " + str(loci_pop2[locus]["a1"]))
        elif loci_pop2[locus]["a1"] >= 10:
            out_f.write("   " + str(loci_pop2[locus]["a1"]))
        
        out_f.write("  " + str(loci_pop2[locus]["a2"]))

print "Done !"
