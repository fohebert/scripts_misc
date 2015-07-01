#!/usr/bin/python

"""

                                \033[1mUser Commmands\033[0m

\033[1mSYNOPSIS\033[0m:
        Takes a fastq file containing reads that have been merged together
        (format for VELVET assembler). The reverse read always follow the
        forward read below. The format is : read 1 FORWARD, then the next
        read in the fastq file is read 1 REVERSE and then same configuration
        for the read #2.
        
        Read_1_F
        Read_1_R
        Read_2_F
        Read_2_R
        ...
        
        This program takes the file and produces 2 output files, one with the
        forward reads and the other file with the reverse reads. 2 files per
        sample (initial fastq input file).

\033[1mUSAGE\033[0m:
        %program <merged_fastq> <forward_file.sfastq> <reverse_file.sfastq>

\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys
import re

try:
    input_file = sys.argv[1]
    output_F = sys.argv[2]
    output_R = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

F = False
R = False
start = 0

with open(input_file, "r") as in_f:
    with open(output_F, "w") as out_F:
        with open(output_R, "w") as out_R:
            for line in in_f:
                start += 1
                line = line.strip()
                if line != "":
                    if line.startswith("@C-"):
                        if re.findall("_F_", line):
                            F = True
                            R = False
                            if start == 1:
                                out_F.write(line)
                            if start > 1:
                                out_F.write("\n" + line)
                    if line.startswith("@C-"):
                        if re.findall("_R_", line):
                            F = False
                            R = True
                            if start == 5:
                                out_R.write(line)
                            if start > 5:
                                out_R.write("\n" + line)
                    else:
                        if F == False and R == True:
                            out_R.write("\n" + line)
                        if F == True and R == False:
                            out_F.write("\n" + line)
