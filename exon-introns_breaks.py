#!/usr/bin/env python

"""

                                    User Commands

\033[1mSYNOPSIS\033[0m:
        Takes a SAM file from a mapping of sequence capture targets on contigs from
        a de novo assembly and outputs the regions of the contigs on which the targets
        mapp. It produces a text file containing, for each contig, the intervals of bp
        for which there is a potential exon.

\033[1mUSAGE\033[0m:
        %program <sam_file_in> <cov_threshold> <exon_length_threshold> <out_file>

\033[1mOUTPUT\033[0m:
        Contig          Exons
        contig_2234     4-50;224-600
        contig_56145    760-880;1000-1089
        ...

\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys
import re

try:
    sam_in = open(sys.argv[1], "rU")
    cov_threshold = int(sys.argv[2])
    length_threshold = int(sys.argv[3])
    output = sys.argv[4]
except:
    print __doc__
    sys.exit(0)

sequences = {}

for line in sam_in:
    
    if line.startswith("@HD"):
        pass
    
    elif line.startswith("@SQ"):
        
        contig = line.split()[1][3:].split("_")[0] + "_" + str("%05i" % int(line.split()[1][3:].split("_")[1]))
        contig_lgth = int(line.split()[2][3:]) + 1
        sequences[contig] = {}
        sequences[contig]["exons"] = {}
        
        for position in range(contig_lgth)[1:]:
            sequences[contig]["exons"][position] = 0
    
    elif line.startswith("@PG"):
        pass
    
    elif line.split()[5].find("S") != -1:
    
        current_contig = line.split()[2].split("_")[0] + "_" + str("%05i" % int(line.split()[2].split("_")[1]))
        
        if re.findall("[0-9]*S[0-9]*M[0-9]*S", line.split()[5]) != []:
            mapped = int(re.findall("[0-9]*M", line.split()[5])[0].split("M")[0])
        
        if re.findall("[0-9]*S[0-9]*M", line.split()[5]) != []:
            unmapped_upstream = int(re.findall("[0-9]*S", line.split()[5])[0].split("S")[0])
            mapped = int(re.findall("[0-9]*M", line.split()[5])[0].split("M")[0])
        
        if re.findall("[0-9]*M[0-9]*S", line.split()[5]) != []:
            mapped = int(re.findall("[0-9]*M", line.split()[5])[0].split("M")[0])
            unmapped_downstream = int(re.findall("[0-9]*S", line.split()[5])[0].split("S")[0])
            
        start_pos = int(line.split()[3])
        end_pos = start_pos + mapped
        
        for exon_pos in range(start_pos,end_pos):
            sequences[current_contig]["exons"][exon_pos] += 1
    
    elif line.split()[5].find("M") and len(line.split()[5]) <= 3:
        
        current_contig = line.split()[2].split("_")[0] + "_" + str("%05i" % int(line.split()[2].split("_")[1]))
        
        mapped = int(line.split()[5].split("M")[0])
        
        start_pos = int(line.split()[3])
        end_pos = start_pos + mapped
        
        for exon_pos in range(start_pos,end_pos):
            sequences[current_contig]["exons"][exon_pos] += 1

with open(output, "w") as out_f:
    
    exon_range = []
    out_f.write("Contig" + "\t" + "Exons")
    
    for contig in list(sorted(sequences)):
        out_f.write("\n" + contig + "\t")
        
        for exon_pos in sequences[contig]["exons"]:
            
            if sequences[contig]["exons"][exon_pos] >= cov_threshold:
                exon_range.append(exon_pos)
            
            if sequences[contig]["exons"][exon_pos] < cov_threshold:                
                
                if exon_range != []:
                    first_pos = sorted(exon_range)[0]
                    last_pos = sorted(exon_range)[len(exon_range) - 1]
                    out_f.write(str(first_pos) + "-" + str(last_pos) + ";")
                    exon_range = []
                else:
                    pass
    
    print "\n\033[1mFile successfully processed\033[0m\n"
