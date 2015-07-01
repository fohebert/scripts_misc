#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
                                  User Commands

\033[1mSYNOPSIS\033[0m

        Formats the file from the script "VCF_parser.py" to create an input file
        for the program BAYESCAN.

\033[1mUSAGE\033[0m

        %program <SNPs_genotypes> <sample_num> <percentage_of_hetero> <output_file>
        
        NB : Percentage of heterozygous allowed to keep for the genotype count. It
        is a FLOAT number between 0.0 and 1.0
        
\033[1mOUTPUT FILE - FORMAT\033[0m
        
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
    in_file = sys.argv[1]
    num_samp = int(sys.argv[2])
    hetero_percent = float(sys.argv[3])
    out_file = sys.argv[4]
except:
    print __doc__
    sys.exit(2)

loci_pop1 = {}
loci_pop2 = {}
smp_col = {}
loci_num = 0
with open(in_file, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        # Associates the name of the sample to its corresponding column position in a dictionnary
        if line.startswith("Contig"):
            samples = num_samp
            ind_column = 5
            
            while samples > 0:
                sp = line.split("\t")[ind_column]
                smp_col[ind_column] = sp
                ind_column += 1
                samples -= 1
        
        elif line.startswith("contig_"):
            
            count = num_samp
            col_num = 5
            pop1_hetero = 0
            pop2_hetero = 0
            while count > 0:
                if line.split("\t")[col_num] != "0":
                    samp_a1 = line.split("\t")[col_num].split("/")[0]
                    samp_a2 = line.split("\t")[col_num].split("/")[1]
                    if samp_a1 == samp_a2:
                        if smp_col[col_num].startswith("C-D"):
                            pop1_hetero += 1
                        if smp_col[col_num].startswith("C-N"):
                            pop2_hetero += 1
                count -= 1
            
            percent_pop1 = float(pop1_hetero)/float(num_samp)
            percent_pop2 = float(pop2_hetero)/float(num_samp)
            
            if percent_pop1 <= hetero_percent and percent_pop2 <= hetero_percent:
                
                loci_num += 1
                loci_pop1[loci_num] = {}
                loci_pop2[loci_num] = {}
                
                tot_num_genes_pop1 = int(line.split("\t")[3].split("/")[0]) + int(line.split("\t")[3].split("/")[1])
                tot_num_genes_pop2 = int(line.split("\t")[4].split("/")[0]) + int(line.split("\t")[4].split("/")[1])
                
                loci_pop1[loci_num]["num_genes"] = 0
                loci_pop1[loci_num]["num_genes"] += tot_num_genes_pop1
                loci_pop2[loci_num]["num_genes"] = 0
                loci_pop2[loci_num]["num_genes"] += tot_num_genes_pop2
                
                loci_pop1[loci_num]["a1"] = 0
                loci_pop1[loci_num]["a1"] += int(line.split("\t")[3].split("/")[0])
                loci_pop1[loci_num]["a2"] = 0
                loci_pop1[loci_num]["a2"] += int(line.split("\t")[3].split("/")[1])
                
                loci_pop2[loci_num]["a1"] = 0
                loci_pop2[loci_num]["a1"] += int(line.split("\t")[4].split("/")[0])
                loci_pop2[loci_num]["a2"] = 0
                loci_pop2[loci_num]["a2"] += int(line.split("\t")[4].split("/")[1])

with open(out_file, "w") as out_f:
    out_f.write("[loci]=" + str(len(loci_pop1)) + "\n" + "\n" + "[populations]=2" + "\n" + "\n" + "[pop]=1")
    
    for locus in sorted(loci_pop1):
        out_f.write("\n")
        
        if locus < 10:
            out_f.write("  " + str(locus) + "  " + str(loci_pop1[locus]["num_genes"]) + "  " + "2")
        elif 10 <= locus < 100:
            out_f.write(" " + str(locus) + "  " + str(loci_pop1[locus]["num_genes"]) + "  " + "2")
        elif locus >= 100:
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
            out_f.write("  " + str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        elif 10 <= locus < 100:
            out_f.write(" " + str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        elif locus >= 100:
            out_f.write(str(locus) + "  " + str(loci_pop2[locus]["num_genes"]) + "  " + "2")
        
        if loci_pop2[locus]["a1"] < 10:
            out_f.write("    " + str(loci_pop2[locus]["a1"]))
        elif loci_pop2[locus]["a1"] >= 10:
            out_f.write("   " + str(loci_pop2[locus]["a1"]))
        
        out_f.write("  " + str(loci_pop2[locus]["a2"]))
