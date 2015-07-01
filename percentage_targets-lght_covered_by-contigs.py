#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                            User Commands

\033[1mSYNOPSIS\033[0m
        This program uses the output file from BLAST+ (format 0) and
        finds the percentage of the length of each target on an exon
        capture  chip that  is covered by the contigs assembled with 
        the chip.

\033[1mUSAGE\033[0m
        %program <blastplus2id> <targets_gene_names>  <blast+> <out_
        contig_IDs> <out>

\033[1mINPUT FILES - FORMAT\033[0m
        
        blastplus2id
        ************
        contig_4	007_006_B1_679590	
        contig_6	CONTIG_MIXED_1043	
        contig_9	CONTIG_MIXED_1185	
        ...
        
        2 columns (tab  delimited) :  contig + corresponding target(s)
        
        This output is produced by the Python script "blastplus2id.py"
        There can be more than one match per contig if the  format is
        respected: the additional matches have to be tab delimited and
        added in different columns to the same target.
        
        \033[1mIt has to be modified though !!!\033[0m
        --> The columns with the e-value must be discarded in order to
            keep only the columns with the target name. Then, a contig
            can have multiple targets, as long as the  format  is  re-
            spected : contig + target1 + target2 + target 3 +  ...
        
        The programs  work  that  way  because it is possible that many
        targets are the same gene, but a different position on the gene
        so when a contig is  assembled, different parts of  the  contig 
        will blast for different targets, which are the same gene.
        
        tagets_gene_names
        *****************
        target1     IGF-B
        target2     S. salar interleukine 1
        target3     heat shock cognate protein 70kD
        ...
        
        blast+
        *******************
        BLAST+ output file, format 0 obtained by blasting the contigs
        against a database of the chip targets.
        
\033[1mOUTPUT FILES - FORMAT\033[0m

        out
        *******************
        Target ID   Gene name     Ref lgth     Contigs     % of ref. lgth
        target1     IGF-B           4587       c1;c24;         45.5
        target2     S. salar...     1495       c324;           79.8
        ...
        
        out_contig_IDs
        *******************
        Contig      Gene name
        contig_1    IGF-B
        contig_2    carboxyesterase
        ...
        
\033[1mCREDIT\033[0m
        Doctor Pants \m/ 2012
"""

import sys
import re

try:
    blastplus2id = sys.argv[1]
    gene_names = sys.argv[2]
    blast_plus_out = sys.argv[3]
    out_contig_IDs = sys.argv[4]
    output_file = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

# Creates a dictionnary with key = contig and value = targets for which
# this contig blasts the best. Will be  used to assign each  contig  to
# the right target when parsing the blast result from BLAST+
contigs = {}
with open(blastplus2id, "rU") as bp_2id:
    for line in bp_2id:
        line = line.strip()
        contig = line.split("\t")[0]
        
        if len(line.split("\t")) > 2:
            num = len(line.split("\t"))
            contigs[contig] = line.split("\t")[1:num]
            
        elif len(line.split("\t")) == 2:
            contigs[contig] = [line.split("\t")[1]]

# Creates a dictionnary  with key = target & value =  gene name for that
# target. Will be used to assign a gene name to each target in the final
# output file.
target_IDs = {}
with open(gene_names, "rU") as g_names:
    for line in g_names:
        line = line.strip()
        target = line.split("\t")[0]
        gene_name = line.split("\t")[1]
        target_IDs[target] = gene_name

targets = {}
contig = ""
current_target = ""
pos_covered = []
good_target = False

# Parses the blast result file from BLAST+ and builds the final dictionnary
# that will be used for the output file.
with open(blast_plus_out, "rU") as blast_f:
    for line in blast_f:
        line = line.strip()
        
        if line.startswith("Query= "):
            
            contig = line.split("Query= ")[1]
            
        if line.startswith("> "):
            
            current_target = line.split("> ")[1]
            
            if current_target in contigs[contig]:
                
                good_target = True
                
                if current_target not in targets:
                    targets[current_target] = {}
                    targets[current_target]["contigs"] = []
                    targets[current_target]["contigs"].append(contig)
                    targets[current_target]["gene_name"] = ""                    
                    targets[current_target]["gene_name"] = target_IDs[current_target]
                    targets[current_target]["lgth_ref"] = 0
                    targets[current_target]["pos_covered"] = set()
                    
                elif current_target in targets:
                    targets[current_target]["contigs"].append(contig)
            
            elif current_target not in contigs[contig]:
                print contig, contigs[contig], current_target
        
        if line.startswith("Length=") and good_target == True:
            
            if targets[current_target]["lgth_ref"] == 0:
                lgth = int(line.split("Length=")[1])
                targets[current_target]["lgth_ref"] += lgth
        
        if re.findall("Sbjct  [0-9]*", line) != [] and good_target == True:
            
            pos = []
            pos.append(int(re.findall("[0-9]* ", line)[2]))
            pos.append(int(re.findall(" [0-9]*", line)[-1]))
            
            for position in range(sorted(pos)[0], sorted(pos)[1]+1):
                targets[current_target]["pos_covered"].add(position)

# Creates the output files
with open(output_file, "w") as out_f:

    out_f.write("Target" + "\t" + "Gene name" + "\t" + "Ref. length" + "\t" + "Contigs" \
        + "\t" + "% of ref. length")
        
    for target in list(sorted(targets)):

        out_f.write("\n" + target + "\t" + targets[target]["gene_name"])
            
        out_f.write("\t" + str(targets[target]["lgth_ref"]) + "\t")
        
        out_f.write(";".join(targets[target]["contigs"]))
        
        percentage = float(len(targets[target]["pos_covered"]))/float(targets[target]["lgth_ref"])
        
        if percentage >= 1.0:
            out_f.write("\t" + "1")
        elif percentage < 1.0:
            out_f.write("\t" + str(percentage))

with open(out_contig_IDs, "w") as out_IDs:
    
    out_IDs.write("Contig" + "\t" + "Gene name")
    
    for contig in list(sorted(contigs)):
        for target in contigs[contig]:
            gene_name = target_IDs[target]
            out_IDs.write("\n" + contig + "\t" + gene_name)

    print "\n\033[1mDone !\033[0m Output file successfully created\n"
