#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
                                        User Commands

\033[1mNAME\033[0m
        dnasp_input - Creates a tab delimited text for import in dnaSP (with SNP alleles)

\033[1mSYNOPSIS\033[0m
        \033[1mdnasp_input.py\033[0m [FILE]...

\033[1mDESCRIPTION\033[0m
        Takes the output SNP table from CLC parsed for the significant SNPs (the ones to
        keep) and based a tab text file containing the contigs associated with the SNPs,
        it creates a FASTA file with 2 sequences for each contig, each sequence contains
        one of the 2 alleles for each SNPs. It can then be imported in dnaSP for analyses
        such as genetic diversity, pN/pS, etc.
        
\033[1mUSAGE\033[0m
        %program <tab_file> <snp_positions> <result_file>
        
        tab_file
        **********************************
        contig00045	CGACCAGAGAATCCATTT...
        contig00666	CAGCTATGTGGTGGATGA...
        contig02182	TTGTGTATTCTTATGACT...
        ...
        
        NOTE : it's possible to convert a FASTA file into a TAB text file with the script
               "\033[1mfasta_to_tab.py\033[0m"
        
        snp_positions
        **********************************
        contig00045  0011  11  SNP  1  A  2  A/T  71,0/29,0  ...
        contig00045  0019  19  SNP  1  A  2  A/T  75,0/25,0  ...
        contig00045  0021  21  SNP  1  G  2  G/A  73,0/27,0  ...
        ...
        
        --> \033[1mWARNING !!!\033[0m <-- : the second column has to have "0"s in front in order to be able
                              to numerically sort the SNP positions properly (the exact same
                              number of numbers for every SNP position).
        
        result_file
        **********************************
        >contig00045_A1
        CGACCAGAGAATCCATTTACGCCATCCAGTACAAC
        >contig00045_A2
        CGACCAGAGATTCCATTTTCACCATCCAGTACAAC
        ...

\033[1mCREDITS\033[0m
        Doctor Pants 2011 \m/
        
"""

import sys

try:
    tab_file = sys.argv[1]
    snp_positions = sys.argv[2]
    result_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

all_contigs = {}
count = 0

with open(snp_positions) as f:
    for line in f:
        line = line.strip()
        if line != "":
            
            count +=1
            
            if count == 1:
                treatment = line.split("\t")[0]
                contig = line.split("\t")[0]
                
                all_contigs[treatment] = {}
                all_contigs[treatment][line.split("\t")[1]] = [line.split("\t")[7][0],
                                                              line.split("\t")[7][2]]
            
            if count > 1:
                contig = line.split("\t")[0]
                
                if contig == treatment:
                    all_contigs[treatment][line.split("\t")[1]] = [line.split("\t")[7][0],
                                                                  line.split("\t")[7][2]]
                
                if contig != treatment:
                    treatment = line.split("\t")[0]
                    all_contigs[treatment] = {}
                    all_contigs[treatment][line.split("\t")[1]] = [line.split("\t")[7][0],
                                                                  line.split("\t")[7][2]]


with open(tab_file) as tab:
    with open(result_file, "w") as out_f:
        for line in tab:
            line = line.strip()
            
            if line != "":
                name = line.split("\t")[0]
                seq = line.split("\t")[1]
                
                sequence1 = seq
                sequence2 = seq
                
                start = 0
                for locus in list(sorted(all_contigs[name])):
                    start += 1
                    
                    if start == 1:
                        
                        position = locus
                        allele1 = all_contigs[name][position][0]
                        allele2 = all_contigs[name][position][1]
                        
                        sequence1 = seq[0:int(position)-1] + allele1
                        sequence2 = seq[0:int(position)-1] + allele2
                        previous_position = locus
                        
                    if start > 1:
                    
                        position = locus
                        allele1 = all_contigs[name][position][0]
                        allele2 = all_contigs[name][position][1]
                        
                        sequence1 += seq[int(previous_position):int(position)-1] + allele1
                        sequence2 += seq[int(previous_position):int(position)-1] + allele2
                        previous_position = locus
                    
                out_f.write(">" + name + "_A1" + "\n" + sequence1 + seq[int(previous_position):] + "\n")
                out_f.write(">" + name + "_A2" + "\n" + sequence2 + seq[int(previous_position):] + "\n")
                
                
                
