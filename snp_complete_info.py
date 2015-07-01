#!/usr/bin/env python

"""
                         User Commands

\033[1mDESCRIPTION\033[0m
\tThis scripts concatenates all the information on SNPs found
\twithin the assembled genes and creates an output  file that
\tcan be used as supplementary information.

\033[1mUSAGE\033[0m
\t%rogram  <SNP_info> <coding_SNPs> <gene_synonymy> <Fst_values> <output>

\033[1mINPUT FILES - DESCRIPTION\033[0m
\tSNP_info --> output from the unifying script
\tgene_synonymy --> info if SNPs are syn or non-syn.

\033[1mCREDITS\033[0m
\tDoctor Pants 2012 \m/
"""

import sys
import re

try:
    SNP_info = sys.argv[1]
    coding_SNPs = sys.argv[2]
    synonymy = sys.argv[3]
    fst_file = sys.argv[4]
    output = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

# Builds up a dictionnary containing the info on coding SNPs ONLY
coding_loci = {}
with open(coding_SNPs, "rU") as coding:
    for line in coding:
        
        gene = line.split("_pos")[0]
        pos = int(re.split("_A[0-9]", line)[0].split("pos")[-1])
        
        # Uses the allele names of the longest allele ORFs and stores the information
        # in the appropriate dictionnary
        if gene not in list(sorted(coding_loci)):
            coding_loci[gene] = set()
            coding_loci[gene].add(pos)
        
        if gene in list(sorted(coding_loci)):
            if pos not in sorted(coding_loci[gene]):
                coding_loci[gene].add(pos)

fsts = {}
snps = {}
with open(synonymy, "rU") as syn_f:
    with open(fst_file, "rU") as fst_f:
        for line in syn_f:
            line = line.strip()
            
            # Establishes the genotype and the position of the coding SNP
            geno = ""
            
            if line.split("\t")[1] != "SNP not included in ORF":
            
                if line.split("\t")[2].split("/")[0][0] != line.split("\t")[2].split("/")[1][0]:
                    geno += line.split("\t")[2].split("/")[0][0] + "" + line.split("\t")[2].split("/")[1][0]
                elif line.split("\t")[2].split("/")[0][1] != line.split("\t")[2].split("/")[1][1]:
                    geno += line.split("\t")[2].split("/")[0][1] + "" + line.split("\t")[2].split("/")[1][1]
                elif line.split("\t")[2].split("/")[0][2] != line.split("\t")[2].split("/")[1][2]:
                    geno += line.split("\t")[2].split("/")[0][2] + "" + line.split("\t")[2].split("/")[1][2]
                
                # Builds the dictionnary that contains the information on the synonymy for
                # coding SNPs
                name = line.split("\t")[0]
                pos = int(line.split("\t")[1])
                
                if name not in list(sorted(snps)):
                    snps[name] = {}
                    snps[name][pos] = {}
                    snps[name][pos]["geno"] = geno
                    snps[name][pos]["synonymy"] = line.split("\t")[-1]
                elif name in list(sorted(snps)):
                    snps[name][pos] = {}
                    snps[name][pos]["geno"] = geno
                    snps[name][pos]["synonymy"] = line.split("\t")[-1]
        
        # Builds the dictionnary that contains the information on Fst values for
        # ALL the SNPs in the dataset
        for line in fst_f:
            line = line.strip()
            
            if line.startswith("Locus_"):
                locus_num = int(line.split("\t")[0].split("Locus_")[-1])
                fst_val = line.split("\t")[2]
                fsts[locus_num] = fst_val

coding = False
with open(SNP_info, "rU") as snp_i:
    with open(output, "w") as out_f:
        
        out_f.write("Gene" + "\t" + "SNP num." + "\t" + "SNP pos." + "\t" + "Coding status" + "\t" + "Geno." \
            + "\t" + "Synonymy" + "\t" + "freq_A1-D" + "\t" + "freq_A1-N" + "\t" + "q-value" + "\t" + "Fst")
        
        for line in snp_i:
            line = line.strip()
            
            if line.startswith("Gene"):
                pass
            
            elif line.startswith("Gene") == False:
                
                coding = False
                name = line.split("\t")[0]
                pos = int(line.split("\t")[1])
                snp_num = int(line.split("\t")[2])
                geno = line.split("\t")[3]
                freq_a1_d = line.split("\t")[4]
                freq_a1_n = line.split("\t")[5]
                all_freq_div = line.split("\t")[6]
                q_value = line.split("\t")[7]
                fst = fsts[snp_num]
                
                out_f.write("\n" + name + "\t" + str(snp_num) + "\t" + str(pos))
                
                if name in list(sorted(coding_loci)):
                    
                    if pos in sorted(coding_loci[name]):
                        out_f.write("\t" + "coding" + "\t" + geno + "\t")
                    elif pos not in sorted(coding_loci[name]):
                        out_f.write("\t" + "non-coding" + "\t" + geno + "\t")
                    
                    if name in list(sorted(snps)):
                        if pos in sorted(snps[name]):
                            out_f.write(snps[name][pos]["synonymy"])
                        elif pos not in sorted(snps[name]):
                            out_f.write("n/a")
                    
                    elif name not in list(sorted(snps)):
                        out_f.write("n/a")
                
                elif name not in list(sorted(coding_loci)):
                    out_f.write("\t" + "non-coding" + "\t" + geno + "\t" + "n/a")
                
                out_f.write("\t" + freq_a1_d + "\t" + freq_a1_n + "\t" + q_value + "\t" + fst)

print "\n\033[1mDone\033[0m\n"
