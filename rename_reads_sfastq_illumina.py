#!/usr/bin/python

"""
Synopsis:
    Rename paired end reads reads with the correct read number and sample name

Usage:
    rename_reads_sfastq_illumina.py  sequence_file  sample_name  output_file

"""

import sys
import re

try:
    in_file = sys.argv[1]
    file_stub = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

read_count_F = 1
read_count_R = 1


with open(in_file, "r") as in_f:
    with open(out_file, "w") as out_f:
        for line in in_f:
            line = line.strip()
            if line != "":
                if line.startswith("@HWI-"):
                    
                    if re.findall("1:N:0:[A-Z]*", line):
                        line = "@C-" + file_stub + "_F_" + str("%07i" % read_count_F)
                        
                        if read_count_F > 1:
                            out_f.write("\n" + line)
                        if read_count_F == 1: 
                            out_f.write(line)
                        
                        read_count_F += 1
                    
                    if re.findall("2:N:0:[A-Z]*", line):
                        line = "@C-" + file_stub + "_R_" + str("%07i" % read_count_R)
                        out_f.write("\n" + line)
                        
                        read_count_R += 1
                
                else:
                    out_f.write("\n" + line)
