#!/usr/bin/env python

"""
                            User Commands

\033[1mDESCRIPTION\033[0m
        This program  takes the output file from the "unifying script"
        and creates, for each exon in which there is a SNP or multiple
        SNPs, a  FASTA  sequence that will be used to find the longest
        ORF  in the  next step.  These ORFs, one for each allele, will
        be used with General  Pick's  script in order to find the syn-
        nymy.

\033[1mUSAGE\033[0m
        %program <all_info_file> <output>

\033[1mOUTPUT FILE - DESCRIPTION\033[0m
        Classic  FASTA file with the name of the gene and the position
        info. (the  nucleotide  position on the gene) + the allele (A1 
        or A2). Example:
        
        >HSP-90_pos234_A1
        ATGGTACCTATGTTATGGCTTATGTTGGGCCTATCTTGATCT
        >HSP-90_pos234_A2
        ATGGTACCTATGTTATGGCTTAAGTTGGGCCTATCTTGATCT
        ...

\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys

try:
    in_file = sys.argv[1]
    out_file = sys.argv[2]
except:
    print __doc__
    sys.exit(2)

exons = {}
with open(in_file, "rU") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            line = line.strip()
            
            count = 1
            if line.split("\t")[3] != "":
                
                gene_name = line.split("\t")[0]
                gene_seq = line.split("\t")[-1]
                all_snp_geno = line.split("\t")[4].split(";")[:-1]
                
                exons = {}
                for exon in line.split("\t")[1].split(";")[:-1]:
                    exons[count] = set()
                    pos1 = int(exon.split("-")[0])
                    pos2 = int(exon.split("-")[-1])
                    
                    for pos in range(pos1, pos2 + 1):
                        exons[count].add(pos)
                    
                    count += 1
                
                snp_count = 0
                for snp_pos in line.split("\t")[3].split(";")[:-1]:
                    
                    snp_pos = int(snp_pos)
                    genotype = all_snp_geno[snp_count]
                    
                    for exon_num in sorted(exons):
                        if snp_pos in sorted(exons[exon_num]):
                            
                            seq_a1 = ""
                            seq_a2 = ""
                            a1 = genotype.split("/")[0]
                            a2 = genotype.split("/")[-1]
                            first_pos = sorted(exons[exon_num])[0]
                            last_pos = sorted(exons[exon_num])[-1]
                            
                            if len(gene_seq) != last_pos:
                                seq_a1 += gene_seq[first_pos - 1:snp_pos - 1] + a1 + gene_seq[snp_pos:last_pos]
                                seq_a2 += gene_seq[first_pos - 1:snp_pos - 1] + a2 + gene_seq[snp_pos:last_pos]
                            
                            elif len(gene_seq) == last_pos:
                                seq_a1 += gene_seq[first_pos - 1:snp_pos - 1] + a1 + gene_seq[snp_pos:]
                                seq_a2 += gene_seq[first_pos - 1:snp_pos - 1] + a2 + gene_seq[snp_pos:]
                            
                            out_f.write(">" + gene_name + "_pos" + str(snp_pos) + "_A1" + "\n" + seq_a1)
                            out_f.write("\n" + ">" + gene_name + "_pos" + str(snp_pos) + "_A2" + "\n" + seq_a2 + "\n")
                    
                    snp_count += 1
