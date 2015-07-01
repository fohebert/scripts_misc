#!/usr/bin/python

"""
v.0.1             User's Command Lines            v.0.1

\033[1mDESCRIPTION\033[0m
    Based on  blast  result files  (concatenated), this
    program   creates a   protein  similarity   profile
    among various species. Proteins for a given species
    was  blasted  (blastp)  against  several  different
    species and the program returns which query protein
    blasted on which species,  according to  thresholds
    as established by Ludin et al. (2011).

\033[1mUSAGE\033[0m
    %program <in_blast> <tag_file> <table_out> <stats_out> <host_proteins>
    
    \033[1mtag_file\033[0m = file with names of species in each category
    Similar to a FASTA file
    
    ________ File example starts below this line ________
    
    >controls
    a_thaliana
    d_melanogaster
    >parasites
    e_granulosus
    e_multilocularis
    >Host
    g_aculeatus
    
    ______________ File example stops here ______________
    
    \033[1mhost_proteins\033[0m = output file containing sequences that blasted on
    the host only.

\033[1mCREDITS\033[0m
    DOC PANTS 2014 \m/
"""

import sys
from Bio import SeqIO
import re

try:
    blast_file = sys.argv[1]
    tag_file = sys.argv[2]
    table_out = sys.argv[3]
    stats_out = sys.argv[4]
    host_proteins = sys.argv[5]
except:
    print __doc__
    sys.exit(1)

# BLASTp identity thresholds placed in a dict() as a reference for later comparison
# Thresholds empirically pre-defined by Ludin et al. 2011
# Threshold dict() for high similarity k-mers (host)
t_denom_pos = (10,11,12,13,14)
t_num_pos = (10,10,11,12,12)
pos_thresholds = dict(zip(t_denom_pos,t_num_pos))

# Threshold dict() for conserved k-mers (controls + parasites)
t_denom_cons = (5,6,7,8,9,10,11,12,13,14)
t_num_cons = (5,5,5,6,7,7,8,8,9,11)
cons_thresholds = dict(zip(t_denom_cons,t_num_cons))

# Dict() created to store the respective category
# of each species, based on the input file given.
# Will be used to count the number of k-mers that
# successfully blasted against species in either
# of the categories given by the user (in this
# case, 3 categories : parasites, controls, host).
# Key = species, value = category. With the species
# name, it's fast and easy to get its category with
# this dictionnary.
species = {}
category = ""
with open(tag_file, "rU") as tag_f:
    for line in tag_f:
        line = line.strip()
        
        if line.startswith(">"):
            category = ""
            category = line.split(">")[-1]
            
        elif line.startswith(">") == False:
            species[line] = category

proteins = {} # Will contain info on kmers per protein and their blast results.
new_query = False
id_start = False
prot = ""
kmer = ""
hit_id = ""
with open(blast_file, "rU") as in_f:
    for line in in_f:
        line = line.strip()
        
        # When encounters a new query, stores k-mer name and
        # protein name and tells the program that a new query
        # was found.
        if line.startswith("Query=") and new_query == False:
            kmer = line.split("Query= ")[-1]
            l = len(kmer.split("_"))
            prot = "_".join(kmer.split("_")[:l-1])
            new_query = True
        
        # In case the new query has no hits found, restarts
        # everything back to zero, ready to find a new query.
        elif re.findall("No hits found", line) != [] and new_query:
            new_query = False
            kmer = ""
            prot = ""
        
        # When a hit for a given query is found
        elif line.startswith("> ") and new_query and id_start == False:
            
            hit_id = "_".join(line.split("> ")[-1].split("_")[:2]) # Stores hit ID
            
            if hit_id in species:
            
                id_start = True # Tells the program there is a new hit
                kmer_num = int(kmer.split("_")[-1]) # Stores k-mer number for the protein
            
                # Will create a new entry in the dict() with the hit ID.
                # Hit ID = species name. Will be used to store info, for 
                # each protein, as to which k-mer blasted on which species
                # with success. At this stage, the dict() is basically getting
                # created and will be filled in when the identity score is found.
                if prot not in proteins:
                    proteins[prot] = {}
                    proteins[prot][kmer_num] = {}
                    proteins[prot][kmer_num][hit_id] = []
                elif prot in proteins:
                    if kmer_num not in proteins[prot]:
                        proteins[prot][kmer_num] = {}
                        proteins[prot][kmer_num][hit_id] = []
                    elif kmer_num in proteins[prot]:
                        proteins[prot][kmer_num][hit_id] = []
            
            # If hit_id is not among the species of interest
            # it won't be included in the dict()
            elif hit_id not in species:
                hit_id = ""
        
        # When identity score is found, stores the info in the big dict()
        elif line.startswith("Identities = ") and new_query and id_start:
            
            ident_num = int(line.split("Identities = ")[1].split(" (")[0].split("/")[0])
            ident_denom = int(line.split("Identities = ")[1].split(" (")[0].split("/")[-1])
            
            # If this k-mer has no blast identity score for this particular
            # species, the scores are entered in the dict(). When hit_id is
            # not empty, this means the k-mer had another hit in this species
            # and it has already been entered in the dict(). Only the first
            # hit per species, i.e. the best, is kept for each k-mer.
            if proteins[prot][kmer_num][hit_id] == []:
                proteins[prot][kmer_num][hit_id].append((ident_num, ident_denom))
            
            hit_id = ""
            id_start = False
        
        # When the query results are over, resets everything back to zero.
        elif line.startswith("Effective search space used:") and new_query:
            new_query = False
            kmer = ""
            prot = ""

# Categories for the Venn diagram are created here.
# Will contain the number of proteins and k-mers
# that blasted against whatever species they
# contain (parasites, controls, host or any combination
# of these categories).
controls = {}
controls["prot"] = 0
controls["kmers"] = 0
parasites = {}
parasites["prot"] = 0
parasites["kmers"] = 0
host = {}
host["prot"] = 0
host["kmers"] = 0
cont_par = {}
cont_par["prot"] = 0
cont_par["kmers"] = 0
cont_par_host = {}
cont_par_host["prot"] = 0
cont_par_host["kmers"] = 0
par_host = {}
par_host["prot"] = 0
par_host["kmers"] = 0
cont_host = {}
cont_host["prot"] = 0
cont_host["kmers"] = 0
host_only = set()
with open(table_out, "w") as t_out:
    with open(stats_out, "w") as s_out:
        
        # Writes down species names in the output
        # file header (tab delimited).
        for spec_name in sorted(species):
            t_out.write("\t" + spec_name)
        
        # Start parsing the big dict() containing info on all
        # proteins.
        for protein in sorted(proteins):
            t_out.write("\n" + protein)
            
            # will contain the num. of k-mer that blasted against 
            # each species for the protein. Start by adding all
            # species names and assigning them the value 0.
            ident_scores = {}
            for spec in species:
                ident_scores[spec] = 0
            
            for kmer in proteins[protein]:
                for hit_id in proteins[protein][kmer]: # hit_id = species name
                    
                    # If species is not the host species, conserved protein thresholds are used
                    if species[hit_id] != "host": # dict() species used to retrieve category
                        if 5 <= proteins[protein][kmer][hit_id][0][1] < 15 :
                            if proteins[protein][kmer][hit_id][0][0] >= cons_thresholds[proteins[protein][kmer][hit_id][0][1]]:
                                ident_scores[hit_id] += 1
                    
                    # If species is host, then high similarity threshold are used.
                    elif species[hit_id] == "host":
                        if proteins[protein][kmer][hit_id][0][1] >= 10:
                            if proteins[protein][kmer][hit_id][0][0] >= pos_thresholds[proteins[protein][kmer][hit_id][0][1]]:
                                ident_scores[hit_id] += 1
            
            in_cont = False
            in_para = False
            in_host = False
            num_kmers = 0
            # For each species name, verifies if there was
            # some k-mers that successfully blasted on it.
            # Verifies which of the 3 categories are
            # represented in the sucessfull blasts (controls,
            # and/or parasites and/or host)
            for spec_name in sorted(ident_scores):
                if ident_scores[spec_name] > 0:
                    if species[spec_name] == "controls":
                        in_cont = True
                        num_kmers += ident_scores[spec_name]
                    elif species[spec_name] == "parasites":
                        in_para = True
                        num_kmers += ident_scores[spec_name]
                    elif species[spec_name] == "host":
                        in_host = True
                        num_kmers += ident_scores[spec_name]
                
                # Writes the number of k-mers that successfully blasted
                # on the species in the "table output file".
                t_out.write("\t" + str(ident_scores[spec_name]))
                
            # fills in the objects destined to create the "stats file" and
            # the Venn diagram (or some sort of inclusion diagram)
            if in_cont == True and in_para == False and in_host == False:
                controls["prot"] += 1
                controls["kmers"] += num_kmers
            elif in_cont == False and in_para == True and in_host == False:
                parasites["prot"] += 1
                parasites["kmers"] += num_kmers
            elif in_cont == False and in_para == False and in_host == True:
                host["prot"] += 1
                host["kmers"] += num_kmers
                host_only.add(protein)
            elif in_cont == True and in_para == True and in_host == False:
                cont_par["prot"] += 1
                cont_par["kmers"] += num_kmers
            elif in_cont == True and in_para == True and in_host == True:
                cont_par_host["prot"] += 1
                cont_par_host["kmers"] += num_kmers
            elif in_cont == False and in_para == True and in_host == True:
                par_host["prot"] += 1
                par_host["kmers"] += num_kmers
            elif in_cont == True and in_para == False and in_host == True:
                cont_host["prot"] += 1
                cont_host["kmers"] += num_kmers
        
        with open(host_proteins, "w") as h_prot:
            for protein in host_only:
                h_prot.write(protein + "\n")
                
        # Writes down, in "stats out file", the number of k-mers in
        # each Venn diagram category. Will be used later to produce
        # a Venn diagram in R.
        s_out.write("Controls only: " + str(controls["prot"]) + " (" + str(controls["kmers"]) + ")")
        s_out.write("\n" + "Parasites only: " + str(parasites["prot"]) + " (" + str(parasites["kmers"]) + ")")
        s_out.write("\n" + "Host only: " + str(host["prot"]) + " (" + str(host["kmers"]) + ")")
        s_out.write("\n" + "Controls and parasites: " + str(cont_par["prot"]) + " (" + str(cont_par["kmers"]) + ")")
        s_out.write("\n" + "Controls, parasites and host: " + str(cont_par_host["prot"]) + " (" + str(cont_par_host["kmers"]) + ")")
        s_out.write("\n" + "Parasites and host: " + str(par_host["prot"]) + " (" + str(par_host["kmers"]) + ")")
        s_out.write("\n" + "Controls and host: " + str(cont_host["prot"]) + " (" + str(cont_host["kmers"]) + ")")
        
print "\n\033[1mJob done !\033[0m\n"