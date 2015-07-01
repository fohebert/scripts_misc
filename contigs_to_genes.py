#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                                User Commands

\033[1mSYNOPSIS\033[0m
        This program is  designed to  concatenate multiple contigs  from
        a de novo assembly into longer  sequences  corresponding to uni-
        genes sequences. It uses the output from a  BLAST+  file (fmt 0)
        to sort  the contigs and associate them to the correct reference
        sequence.
        
        The BLAST+ file should be the contigs  blasted  against a  local
        database corresponding to the reference sequences (unigenes).

\033[1mUSAGE\033[0m
        %program <FASTA_contigs> <gene_names> <blast+_result> <output>
            <output_no_blast>

\033[1mINPUT FILE - FORMAT\033[0m
        gene_names
        **********
        List of the gene names composing the reference sequences,one name
        per line.

\033[1mOUTPUT FILE - FORMAT\033[0m
        FASTA file containing the concatenated contigs with the names of
        the reference sequences (i.e gene names).
        
        output_no_blast
        ***************
        List of ref sequences (ref genes) that didn't get any blast hit.
        One name per line in a normal text file.

\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

try:
    contigs_fasta = open(sys.argv[1], "rU")
    gene_names = sys.argv[2]
    blast = sys.argv[3]
    output = sys.argv[4]
    out_no_hit = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(contigs_fasta, "fasta"))

# Fills in a dictionnary that contains the contig names with their corresponding sequence
# Will be used  later on to  reverse complement certain sequences and to re-construct the
# gene sequences
contigs = {}
for seq in sequences:
    contig = seq[0]
    contig_seq = seq[1]
    contigs[contig] = {}
    contigs[contig]["rev"] = ""
    contigs[contig]["original"] = ""
    contigs[contig]["original"] += contig_seq

# Starts constructing the dictionnary containing the info on each gene,
# i.e the contigs and their position on the sequence of those ref genes.
genes = {}
with open(gene_names, "rU") as g_names:
    for line in g_names:
        line = line.strip()
        gene = line
        genes[gene] = {}

# Finishes the dictionnary with the info on each ref gene.
temporary = {}
contig = ""
gene = ""
first_query = False
started = False
new_gene = False
with open(blast, "rU") as b_file:
    for line in b_file:
        line = line.strip()
        
        if line.startswith("Query="):
            
            if first_query == False:
                contig = line.split("Query= ")[1]
                first_query = True
                
            elif first_query:
                
                if len(temporary) > 1:
                    
                    # Verifies among the genes if there is a degenerate paralog sequence
                    # If the contig blasts for several genes among which there is a paralog
                    # sequence, it pops it out of the dictionnary.
                    for gene in sorted(temporary):
                        
                        if len(temporary) > 1:
                        
                            if re.findall("probable degenerate paralog", gene):
                                temporary.pop(gene)
                    
                    # If there is still more than one gene for which the contig has blasted,
                    # it keeps the hit that has the higher score and discards the other one.
                    # In the end, only one gene is kept of the contig.
                    last_best_gene = ""
                    blast_score = 0
                    for gene in temporary:
                        if last_best_gene == "":
                            last_best_gene = gene
                            blast_score += temporary[gene]["score"]
                        elif last_best_gene != "":
                            if temporary[gene]["score"] > blast_score:
                                last_best_gene = gene
                                blast_score = temporary[gene]["score"]
                            elif temporary[gene]["score"] <= blast_score:
                                pass
                    
                    # Stores the hit positions of the contig on the ref. gene in the object "pos"
                    
                    pos = str(list(sorted(temporary[last_best_gene]["positions"]))[0]) + "-" \
                        + str(list(sorted(temporary[last_best_gene]["positions"]))[-1])
                    
                    # Fills in the final dictionnary : associates the contig with the right gene
                    # and stores the positions on the gene for which the contig has blasted plus
                    # the blast score. This dictionnary will be used later to build the final
                    # FASTA file.
                    genes[last_best_gene][contig] = {}
                    genes[last_best_gene][contig]["score"] = blast_score
                    genes[last_best_gene][contig]["positions"] = pos
                
                elif len(temporary) == 1:
                    
                    for gene in temporary:
                        pos = str(list(sorted(temporary[gene]["positions"]))[0]) + "-" \
                            + str(list(sorted(temporary[gene]["positions"]))[-1])
                        
                        genes[gene][contig] = {}
                        genes[gene][contig]["score"] = temporary[gene]["score"]
                        genes[gene][contig]["positions"] = pos
                
                contig = line.split("Query= ")[1]
                temporary = {}
        
        # Stores the name of the gene for which the contig blasted in the object "gene"
        elif line.startswith("> "):
            gene = line.split("> ")[1]
            started = False
            new_gene = True
        
        # Fills in the dictionnary "temporary" with the blast info for each gene (ref seq)
        # for which the contig blasts. Sorts and stores the info according to the blast
        # score, which will be used later to determine which "true" gene the contig
        # represents
        elif re.findall("Score =", line) != []:
            
            if new_gene:
            
                # If this is a new gene associated to the contig, it fills in the
                # dictionnary containing the necessary temporary info on that blast
                # result 
                score = int(line.split()[2])
                temporary[gene] = {}
                temporary[gene]["score"] = score
                temporary[gene]["positions"] = set()
            
            # If it is the same gene, but a different position on this gene, it just completes
            # the info in the temporary dictionnary concerning the blast positions on the ref
            # seq.
            elif new_gene == False:
                
                # It adds the score to the total score for that gene in the dictio. "temporary"
                score = int(line.split()[2])
                temporary[gene]["score"] += score
                
            new_gene = False
        
        elif line.startswith("Sbjct"):
            if started:
                # Adds the positions on the reference for which the contig blasts
                # in a dictionnary that will be used later for the final info 
                pos1 = int(re.findall("Sbjct  [0-9]*", line)[0].split()[1])
                pos2 = int(re.findall("[0-9]*$", line)[0])
                temporary[gene]["positions"].add(pos1)
                temporary[gene]["positions"].add(pos2)
            
            elif started == False:
                started = True
                # Adds the positions on the reference for which the contig blasts
                # in a dictionnary that will be used later for the final info 
                pos1 = int(re.findall("Sbjct  [0-9]*", line)[0].split()[1])
                pos2 = int(re.findall("[0-9]*$", line)[0])
                temporary[gene]["positions"].add(pos1)
                temporary[gene]["positions"].add(pos2)
                # Reverse comnplement the contig sequence if its orientation
                # is opposed to that of the reference sequence
                if pos1 > pos2:
                    if contigs[contig]["rev"] == "":
                        seq_ori = contigs[contig]["original"]
                        dna = Seq(seq_ori, generic_dna)
                        seq_rev = dna.reverse_complement().tostring()
                        contigs[contig]["rev"] += seq_rev

# Generates the output file
with open(output, "w") as out_f:
    with open(out_no_hit, "w") as no_hit:
        for gene in list(sorted(genes)):
            
            if genes[gene] == {}:
            
                no_hit.write(gene + "\n")
            
            elif genes[gene] != {}:
            
                out_f.write(">" + gene)
                gene_positions = {}
                to_discard = set()
                contigs_to_compare = []
                
                for contig in genes[gene]:
                    
                    # Uses the dictionnary "gene_positions" to store the information on every contig
                    # blasting for a certain gene. Will be used to discard the contigs that are falsely
                    # associated to that gene and to order the good contigs to build the final whole
                    # sequence.
                    tot_pos = genes[gene][contig]["positions"]
                    gene_positions[contig] = {}
                    gene_positions[contig]["tot_pos"] = tot_pos
                    gene_positions[contig]["score"] = genes[gene][contig]["score"]
                    contigs_to_compare.append(contig)
                
                # Will use the list of all the contigs for a specific gene to compare them among each
                # other. It compares the overlapping positions of every contig with all the other
                # contigs for that gene to see if there are some that should be discarded.
                for current_contig in contigs_to_compare:
                    
                    # If the contig that we are comparing to the others is still flagged as
                    # "to discard", it is not compared. If it hasn't been compared yet or if
                    # it doesn't overlap with any contig yet, it will be compared to the others.
                    if current_contig not in to_discard:
                    
                        # creates a set with the blast positions on the ref gene for the current contig
                        locus = set()
                        current_pos1 = int(gene_positions[current_contig]["tot_pos"].split("-")[0])
                        current_pos2 = int(gene_positions[current_contig]["tot_pos"].split("-")[-1])
                        for base in range(current_pos1, current_pos2 + 1):
                            locus.add(base)
                        
                        # Iterates over all the contigs that blast for the same gene giving that
                        # the contigs are not already flagged as "to_discard".
                        for contig in gene_positions:
                            
                            if current_contig == contig:
                                pass
                            
                            elif current_contig != contig and contig not in to_discard:
                                
                                count = 0
                                contig_pos1 = int(gene_positions[contig]["tot_pos"].split("-")[0])
                                contig_pos2 = int(gene_positions[contig]["tot_pos"].split("-")[-1])
                                for i in range(contig_pos1, contig_pos2 + 1):
                                    if i in locus:
                                        count += 1
                                
                                if count >= 30:
                                    if gene_positions[current_contig]["score"] > gene_positions[contig]["score"]:
                                        to_discard.add(contig)
                                    if gene_positions[current_contig]["score"] < gene_positions[contig]["score"]:
                                        to_discard.add(current_contig)
                                    if gene_positions[current_contig]["score"] == gene_positions[contig]["score"]:
                                        if len(range(contig_pos1, contig_pos2 + 1)) > len(locus):
                                            to_discard.add(current_contig)
                                        else:
                                            to_discard.add(contig)
                
                # After all the iterations and the comparisons between the contigs, the ones that have to be
                # discarded in the final sequence are taken out of the dictionnary "gene_positions"
                for contig_name in to_discard:
                    gene_positions.pop(contig_name)
                
                # The dictionnary contig_order is created to sort the contigs and put them in the correct
                # order to construct the final sequence of the gene with the contigs that have been kept
                # after the last loop and the comparisions among the contigs altogether.
                contig_order = {}
                for contig in gene_positions:
                
                    # Keeps only the first position of the region covered by the contig on the ref gene and
                    # sorts them according to that. The contig with the smallest first position is the first
                    # contig of the final sequence and so on. 
                    first_pos = int(gene_positions[contig]["tot_pos"].split("-")[0])
                    contig_order[first_pos] = contig
                
                # Whole sequence is the final sequence of the gene. By sorting the entries in the dictionnary 
                # "contig_order" the contigs are successively added to the object "whole_sequence".
                whole_sequence = ""
                for first_pos in sorted(contig_order):
                    
                    contig = contig_order[first_pos]
                    
                    if contigs[contig]["rev"] != "":
                        contig_seq = contigs[contig]["rev"]
                        whole_sequence += contig_seq
                    
                    elif contigs[contig]["rev"] == "":
                        contig_seq = contigs[contig]["original"]
                        whole_sequence += contig_seq
                    
                # When the loop is over and all the contig sequences have been added successively, the object
                # "whole_sequence" is written in the output file in the FASTA format (in one line though).
                out_f.write("\n" + whole_sequence + "\n")
        
    print "\nOutput file successfully created\n"
