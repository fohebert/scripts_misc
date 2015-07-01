#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                            User's Commands

\033[1mSYNOPSIS\033[0m
        Takes a .gff file from the program AUGUSTUS and returns a file
        with the predicted gene names, sequence (a.a.),  saffold  name
        and scaffold position (on reference genome). It  also produces
        a FASTA file with the protein sequences of each predicted gene.

\033[1mUSAGE\033[0m
        %program <GFF_input_file> <output_table> <output_FASTA>
        
\033[1mOUTPUT FILE - FORMAT\033[0m

         _________________file starts below the line_________________
                
         Gene                Scaffold        Position      Seq (a.a.)
         pred_gene_0001       00001         34234-34600    KSTNKKN...
         pred_gene_0002       00001         62004-62739    AVGSSST...
                
         ___________________file ends above the line_________________

\033[1mCREDITS\033[0m
        Doctor Pants 2013 \m/
"""

import sys
import re

try:
    in_gff = sys.argv[1]
    out_file = sys.argv[2]
    out_fasta = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

g_prefix = raw_input('\nWhat kind of prefix do you want for your genes ? Enter here --> ')

header = ['Gene', 'Scaffold', 'Position', 'Sequence (a.a.)']
scaffold = ""
seq = ""
start_seq = False
with open(in_gff, "rU") as i_gff:
    with open(out_fasta, "w") as out_fas:
        with open(out_file, "w") as out_f:
            
            out_f.write("\t".join(header))
            
            for line in i_gff:
                line = line.strip()
                
                # Extracts the current scaffold name when a new one is encountered
                if line.startswith('# ----- prediction on sequence number'):
                    scaffold = ""
                    scaffold = line.split("name = ")[-1].split(") -----")[0]
                
                # Extracts the current gene name when a new one is encountered and writes 
                # it down in the output file. All  genes  will  have  the prefix specified 
                # in the command line arguments.
                elif line.startswith('# start gene'):
                    gene = g_prefix + ('%04i' % int(line.split()[-1].split('g')[-1]))
                    out_f.write('\n' + gene + '\t' + scaffold)
                    out_fas.write('>' + gene) # Writes down the gene name in the FASTA output
                
                # Extracts the position of the current gene on the current scaffold and
                # writes it down in the output file.
                elif line.startswith(scaffold):
                    if line.find('\tgene\t') != -1:
                        out_f.write("\t" + line.split("\t")[3] + "-" + line.split("\t")[4])
                
                # Extracts the protein sequence (a.a.) of the current gene (first line)
                elif line.startswith('# protein sequence'):
                    start_seq = True
                    seq += line.split(" [")[-1]
                
                # Extracts subsequent chunks of current protein sequence
                elif line.startswith('# ') and start_seq:
                    # Puts the fragment of the protein sequence in the object 'seq', 
                    # line by line until it reaches the end, marked by ']'
                    if line.find(']') != -1:
                        seq += line.split('# ')[-1].split(']')[0]
                        start_seq = False
                        out_f.write('\t' + seq)
                        
                        # Loop that takes the final sequence contained in the object
                        # 'seq' and writes it down in FASTA format in the appropriate
                        # output file. FASTA file contains sequences sorted by gene
                        # name, 60 ncl per line.
                        for i in range(0,len(seq),60):
                            chunk = seq[i:i+60]
                            out_fas.write('\n' + chunk)
                        
                        out_fas.write('\n')
                        
                        # Empties 'seq' object so that it can be filled with another
                        # sequence (new one)
                        seq = ""
                        
                    if line.find(']') == -1:
                        seq += line.split('# ')[-1]

print "\n\033[1mDone\n"
