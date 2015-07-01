#!/usr/bin/env python

"""
                            User Commands

\033[1mSYNOSPIS\033[0m
        Uses a sam file from a mapping of reads on reference sequences
        and outputs a file that can be used to produce a graph of type
        "cumulative  base frequency".  This  graph will help to deter-
        mine if  the coverage PER BASE is uniform and sufficient among 
        all the reference sequences.

\033[1mUSAGE\033[0m
        %program <sam_file> <tot_num_bases_in_ref_seq> <output> <log>

\033[1mOUTPUT FILE - FORMAT\033[0m
        Coverage    Cumulative base frequency
        1           97.7
        2           97.7
        3           97.7
        ...         ...
        10000       1.2
        
        LOG FILE
        ********
        Contains the information on the dictionnaries within the script
        to verify if the final percentage is correct.
        
\033[1mCREDITS\033[0m
        Doctor Pants 2012 \m/
"""

import sys
import re

try:
    sam_file = sys.argv[1]
    tot_bases = int(sys.argv[2])
    out_file = sys.argv[3]
    log = sys.argv[4]
except:
    print __doc__
    sys.exit(0)

ref_seq = {}
with open(sam_file, "rU") as sam_f:
    for line in sam_f:
        line = line.strip()
            
        if line.startswith("@HD"):
            pass
        
        # Writes in the big dictionnary the name of the reference sequences that are in the SAM header
        elif line.startswith("@SQ"):
            ref_seq_name = line.split("\t")[1].split(":")[1]
            ref_seq[ref_seq_name] = {}
                
        elif line.startswith("@PG"):
            pass
        
        # Takes the info for all the reads, i.e their position on the ref. seq.
        # It enters in the big dictionnary the coverage + position information
        elif line.startswith("HWI"):
            
            # In the case the whole read maps to the ref. seq. (no INDELs)
            if re.findall("100M", line) != []:
                
                ref_seq_name = line.split("\t")[2]
                first_pos = int(line.split("\t")[3])
                last_pos = first_pos + 100
                
                # Goes through every position on the ref. seq. from the first one to the last one    
                for pos in range(first_pos,last_pos):
                    
                    # Adds 1 to the position on the ref seq., indicating that there is a read for
                    # that position (it counts the coverage, position by position)
                    if pos not in ref_seq[ref_seq_name]:
                        ref_seq[ref_seq_name][pos] = 0
                        ref_seq[ref_seq_name][pos] += 1
                    elif pos in ref_seq[ref_seq_name]:
                        ref_seq[ref_seq_name][pos] += 1
            
            # In the case the read does not fully map to the ref. seq. (gaps or INDELs)
            elif re.findall("[0-9]*[A-Z][0-9]*[A-Z]", line.split("\t")[5]) != []:
                
                ref_seq_name = line.split("\t")[2]
                first_pos = int(line.split("\t")[3])
                current_pos = first_pos
                
                # Separates the maping pattern, i.e which positions map, which positions are gaps and
                # which positions are insertions.
                map_pattern = re.findall("[0-9]*[A-Z]", line.split("\t")[5])
                for i in map_pattern: # goes through each pattern
                        
                    if re.split("[0-9]*", i)[1] == "M": # M = map
                        
                        num_to_add = int(re.split("[A-Z]", i)[0])
                        
                        if current_pos == first_pos:
                            for pos in range(current_pos, current_pos + num_to_add):
                                    
                                if pos not in ref_seq[ref_seq_name]:
                                    ref_seq[ref_seq_name][pos] = 0
                                    ref_seq[ref_seq_name][pos] += 1
                                elif pos in ref_seq[ref_seq_name]:
                                    ref_seq[ref_seq_name][pos] += 1
                                
                            current_pos += num_to_add - 1
                            
                        elif current_pos != first_pos:
                            for pos in range(current_pos + 1, current_pos + num_to_add + 1):
                                    
                                if pos not in ref_seq[ref_seq_name]:
                                    ref_seq[ref_seq_name][pos] = 0
                                    ref_seq[ref_seq_name][pos] += 1
                                elif pos in ref_seq[ref_seq_name]:
                                    ref_seq[ref_seq_name][pos] += 1
                                
                            current_pos += num_to_add
                                    
                    elif re.split("[0-9]*", i)[1] == "I": # I = insertion
                        pass
                        
                    elif re.split("[0-9]*", i)[1] == "D": # D = deletion
                        
                        num_to_add = int(re.split("[A-Z]", i)[0])
                        current_pos += num_to_add
                            
                    elif re.split("[0-9]*", i)[1] == "S":
                        pass

num_reads = {}
with open(log, "w") as log:
    for seq in list(sorted(ref_seq)):
        for pos in ref_seq[seq]:
            
            # For every position in every ref. seq. it adds 1 to the corresponding
            # coverage category for every read representing that base. If there is
            # 10 reads, it will add 1 in each coverage category, from 1 to 10 and
            # so one for every position of every ref. seq.
            for i in range(1, ref_seq[seq][pos] + 1):
                
                if i not in num_reads:
                    num_reads[i] = 0
                    num_reads[i] += 1
                
                elif i in num_reads:
                    num_reads[i] += 1

    log.write("Number of targets = " + str(len(ref_seq)) + "\n")

    log.write("Coverage counts" + "\n")
    for cov in num_reads:
        log.write(str(cov) + "\t" + str(num_reads[cov]) + "\n")
    
    log.write("Targets + num_pos_with_reads")
    for seq in list(sorted(ref_seq)):
        log.write(seq + "\t" + str(len(ref_seq[seq])) + "\n")

with open(out_file, "w") as out_f:
    out_f.write("Coverage" + "\t" + "Cumulative base frequency")
    
    # Simply writes down in the output file the number of bases for each coverage category
    # The info was written in the dictionnary "num_reads" filled in in the previous step.
    for cov in num_reads:
        out_f.write("\n" + str(cov) + "\t")
        
        frequency = (float(num_reads[cov])/float(tot_bases))*100
        out_f.write(str(frequency))
