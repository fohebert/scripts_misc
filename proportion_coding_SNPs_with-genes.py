#!/usr/bin/env python

"""
                                    User Commands

\033[1mSYNOPSIS\033[0m
        This program uses the final output of the script "complete_contig_info.py"
        and  returns on  the  screen the number of coding and non-conding SNPs. It
        also  creates  2 output files each containing the contig sequences for one 
        of the 2 SNP alleles. These output files can then be used to calculate dn/
        ds ratios with MEGA. 
        
        Format of the output files = FASTA.
        
\033[1mUSAGE\033[0m
        %program <complete_contig-info> <output_prefix> <dnasp_out> <list_genes>
        
\033[1mCREDITS\033[0m
        Doc Pants 2012 \m/
"""

import sys
import re

try:
    all_info = sys.argv[1]
    out_prefix = sys.argv[2]
    dnasp_out = sys.argv[3]
    list_genes = sys.argv[4]
except:
    print __doc__
    sys.exit(2)

num_coding = 0
num_non_coding = 0
total_snps = 0
gene_count = 1
num_contigs = 0
coding_contigs = set()
coding = False
with open(all_info, "rU") as in_f:
    with open(dnasp_out, "w") as dnasp:
        with open(list_genes, "w") as l_genes:
            for line in in_f:
                line = line.strip()
            
                if line.startswith("Gene"):
                    pass
            
                if line.startswith("Gene") == False:
                    
                    coding = False
                    seq1 = ""
                    seq2 = ""
                    if re.findall("[A-Z]/[A-Z];", line) != []:
                        coding_region = set()
                        contig_name = line.split("\t")[0]
                        num_contigs += 1
                        
                        for exon in line.split("\t")[1].split(";"):
                            if exon != "No exons found" and exon != "":
                                pos1 = int(exon.split("-")[0])
                                pos2 = int(exon.split("-")[1])
                                for i in range(pos1,pos2+1):
                                    coding_region.add(i)
                        
                        snp_pos = {}
                        count = 0
                        for pos in line.split("\t")[3].split(";"):
                            if pos != "":
                                pos = int(pos)
                                if pos in coding_region:
                                    num_coding += 1
                                    total_snps += 1
                                    coding_contigs.add(contig_name)
                                    snp_pos[pos] = line.split("\t")[4].split(";")[count]
                                    coding = True
                                elif pos not in coding_region:
                                    num_non_coding += 1
                                    total_snps += 1
                                count += 1
                        
                        lgth = 0
                        for pos in sorted(coding_region):
                            lgth += 1
                            
                            if pos in sorted(snp_pos):
                                seq1 += snp_pos[pos].split("/")[0]
                                seq2 += snp_pos[pos].split("/")[-1]
                            elif pos not in sorted(snp_pos):
                                seq1 += line.split("\t")[-1][pos - 1]
                                seq2 += line.split("\t")[-1][pos - 1]
                        
                        if seq1 != "" and seq2 != "" and coding:
                            with open(out_prefix + "_gene" + str("%04i.fas" % gene_count), "w") as out_f:
                                out_f.write(">" + contig_name + "_A1" + "\n" + seq1)
                                out_f.write("\n" + ">" + contig_name + "_A2" + "\n" + seq2 + "\n")
                            
                            l_genes.write(contig_name + "\t" + str(len(coding_contigs)) + "\n")
                            
                            dnasp.write(">" + contig_name + "_A1" + "\n" + seq1 + "\n" + ">" + \
                                ">" + contig_name + "_A2" + "\n" + seq2 + "\n")
                            
                            seq1 = ""
                            seq2 = ""
                        
                    gene_count += 1
            
    print "\nNumber of coding SNPs: ", str(num_coding), "(",len(coding_contigs)," genes)\n"
    print "Number of non coding SNPs: ", str(num_non_coding), "\n"
    print "Total number of SNPs: ", str(total_snps), "(",num_contigs," genes)\n"
