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
        %program <complete_contig-info> <contigs-allele1> <contigs-allele2>
        
\033[1mCREDITS\033[0m
        Doc Pants 2012 \m/
"""

import sys
import re

try:
    all_info = sys.argv[1]
    out1 = sys.argv[2]
    out2 = sys.argv[3]
except:
    print __doc__
    sys.exit(2)

num_coding = 0
num_non_coding = 0
total_snps = 0
with open(all_info, "rU") as in_f:
    with open(out1, "w") as out_f1:
        with open(out2, "w") as out_f2:
            
            for line in in_f:
                line = line.strip()
                
                if line.startswith("Contig"):
                    pass
                
                if line.startswith("contig_"):
                    
                    if re.findall("[A-Z]/[A-Z];", line) != []:
                        coding_region = set()
                        contig_name = line.split("\t")[0]
                        out_f1.write(">" + contig_name + "\n")
                        out_f2.write(">" + contig_name + "\n")
                        
                        for exon in line.split("\t")[2].split(";"):
                            if exon != "":
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
                                elif pos not in coding_region:
                                    num_non_coding += 1
                                    total_snps += 1
                                    
                                snp_pos[pos] = line.split("\t")[4].split(";")[count]
                                count += 1
                        
                        seq1 = ""
                        seq2 = ""
                        contig_pos = 0
                        for pos in list(sorted(snp_pos)):
                            seq1 += line.split("\t")[-1][contig_pos:pos-1] + snp_pos[pos].split("/")[0]
                            seq2 += line.split("\t")[-1][contig_pos:pos-1] + snp_pos[pos].split("/")[1]
                            contig_pos = pos
                        
                        if contig_pos < len(line.split("\t")[-1]):
                            seq1 += line.split("\t")[-1][pos:]
                            seq2 += line.split("\t")[-1][pos:]
                        
                        out_f1.write(seq1 + "\n")
                        out_f2.write(seq2 + "\n")
            
            print "\nNumber of coding SNPs: ", str(num_coding), "\n"
            print "Number of non coding SNPs: ", str(num_non_coding), "\n"
            print "Total number of SNPs: ", str(total_snps), "\n"
