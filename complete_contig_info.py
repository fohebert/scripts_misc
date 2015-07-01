#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

                               User Commands

\033[1mSYNOPSIS\033[0m
        This program gathers all the necessary information on the contigs
        captured by an exon capture chip from multiple files and outputs 
        it in a text file in which each contig (one contig per line) is
        associated with various characteristics, i.e gene name, exon po-
        sitions, SNP positions, SNP genotypes and contig sequence.

\033[1mUSAGE\033[0m
        %program <contig_IDs> <fasta_file> <exon_file> <SNP_genotypes> <output_file>

\033[1mINPUT FILES FORMAT\033[0m
        <contig_IDs>
        One contig name per line + tab + gene name, ex. : 
        Contigs        Gene name
        contig_1045    IGF-B
        contig_0021    Heat Shock Cognate Protein 70K
        ...
        
        <fasta_file>
        Usual fasta file for all the contigs analysed
        
        <exon_file>
        Output file from the script "exon-intron_breaks".
        
        <SNP_genotypes>
        Output file from the script "VCF_parser.py"
        
\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/

"""

import sys
from Bio import SeqIO

try:
    gene_name = sys.argv[1]
    fasta_file = open(sys.argv[2], "rU")
    exon_file = sys.argv[3]
    SNP_genotypes = sys.argv[4]
    output_file = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

contigs = {}

sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(fasta_file, 'fasta'))

# Creates the meta-dictionnary and fills out each contig parameter with empty strings
# except for the sequence of the contig, which is added to the dictionnary base on the
# FASAT file given in the arguments of the command line.
for seq in sequences:
    contig_name = seq[0].split("_")[0] + "_" + str("%05i" % int(seq[0].split("_")[1]))
    contigs[contig_name] = {}
    contigs[contig_name]["sequence"] = seq[1]
    
    contigs[contig_name]["exons"] = ""
    contigs[contig_name]["contig_ID"] = ""
    contigs[contig_name]["snps"] = {}

print "Sequences processed\nStarting to assign the information to each contigs..."

# Fills in the contig_ID parameter
gene = {}
with open(gene_name, "rU") as gene_n:
    for line in gene_n:
        line = line.strip()
        
        if line.startswith("Contig"):
            pass
        
        elif line.startswith("contig_"):
            
            contig_name = line.split("\t")[0].split("_")[0] + "_" + str("%05i" % int(line.split("\t")[0].split("_")[1]))
            gene_name = line.split("\t")[1]
            
            if contig_name not in gene:
                gene[contig_name] = []
                gene[contig_name].append(gene_name)
            elif contig_name in gene:
                gene[contig_name].append(gene_name)
    
    for contig in list(sorted(gene)):
        
        contigs[contig]["contig_ID"] = gene[contig]
        
#        if len(line.split("\t")) > 1:
#            contig = line.split("\t")[0].split("_")[0] + "_" + str("%05i" % int(line.split("\t")[0].split("_")[1]))
#            contigs[contig]["contig_ID"] = line.split("\t")[1]
#        if len(line.split("\t")) == 1:
#            contig = line.split("\t")[0]
#            contigs[contig]["contig_ID"] = "unknown gene ID"
        
    print "\nGene IDs found and processed"

# Fills in the exon positions parameter
with open(exon_file, "rU") as exon_f:
    for line in exon_f:
        line = line.strip()
        if line.startswith("Contig\t"):
            pass
        elif line.startswith("contig_"):
            if len(line.split()) > 1:
                contig = line.split()[0]
                exons = line.split()[1]
                contigs[contig]["exons"] = exons
            if len(line.split()) == 1:
                contig = line.split()[0]
                contigs[contig]["exons"] = "-"
    
    print "\nExon position processed..."

# Fills in the SNPs parameters (postion AND genotype)
with open(SNP_genotypes, "rU") as snp_f:
    for line in snp_f:
        line = line.strip()
        if line.startswith("Contig\t"):
            pass
        elif line.startswith("contig_"):
            contig = line.split()[0]
            contigs[contig]["snps"][int(line.split()[1])] = line.split()[2]

    print "\nSNPs processed..."

# Creates the output file
with open(output_file, "w") as out_f:
    
    # Creates the file header
    out_f.write("Contig" + "\t" + "Gene ID" + "\t" + "Exons" + "\t" \
        "SNP position" + "\t" + "SNP genotype" + "\t" + "Contig sequence")
    
    for contig in list(sorted(contigs)):
        out_f.write("\n" + contig + "\t")
        
        out_f.write("---".join(contigs[contig]["contig_ID"]) + "\t" + contigs[contig]["exons"] + "\t")
        
#        out_f.write(contigs[contig]["contig_ID"] + "\t" + contigs[contig]["exons"] + "\t")
        
        for snp in list(sorted(contigs[contig]["snps"])):
            out_f.write(str(snp) + ";")
        
        out_f.write("\t")
        
        for snp in list(sorted(contigs[contig]["snps"])):
            out_f.write(contigs[contig]["snps"][snp] + ";")
        
        out_f.write("\t" + contigs[contig]["sequence"])
    
    print "\nDONE - Output file successfully created\n"
