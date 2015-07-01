#!/usr/bin/python

"""
\033[1mDESCRIPTION\033[0m
    Coverage calculator. Will take a list of contig name with their number of reads
    (column 1 = contig name & column 2 = read count) and will calculate, for each
    contig, the coverage. This program needs an important prior: the average read
    length.

    Output file will look like this: one contig per line, tab-delimited, with first
    column = contig name and second column = coverage for that contig. The program
    will print on the screen the average coverage for the whole assembly, i.e. the
    average for all of these individual values written in the output file.

\033[1mUSAGE\033[0m
    %program <list_read_count> <fasta> <output> <average read length>

\033[1mCREDITS\033[0m
    Doc Pants 2015 \m/
"""

import sys
from Bio import SeqIO

try:
    read_count = sys.argv[1]
    fasta_in = open(sys.argv[2], "rU")
    output = sys.argv[3]
    ave_read_lgth = float(sys.argv[4])
    
except:
    print __doc__
    sys.exit(215)

##################
# GLOBAL VARIABLES
##################

# Creates an object called "sequences" that contain the name + the sequence of each contig
# in the assembly. The info can be accessed through seq[0] = seq name & seq[1] = nucl. seq.
sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(fasta_in, "fasta"))

# Dictionnary containing all of the contig lengths
lengths = {}
for seq in sequences:
    lengths[seq[0]] = float(len(seq[1]))

# List that will contain all of the coverages calculated and will be used to calculate
# the average coverage for the whole assembly.
coverages = []

###########
# FUNCTIONS
###########

# Finds the median in a list (i.e. list containing all of the coverage values)
def median(lst):
    half = len(lst)/2
    lst.sort()
    if len(lst) % 2 == 0: # if list length is even, average middle two values
        return (lst[half-1] + lst[half])/2
    else:
        return lst[half]

# Compute the range of values 
def range(lst):
    lst.sort()
    return [lst[0],lst[-1]]

# Coverage computation
def cov(r_count,ave_r_lgth,contig_lgth):
    return (float(r_count)*float(ave_r_lgth)) / float(contig_lgth)

# Average coverage computation
def ave_cov(lst):
    return float(sum(lst))/float(len(lst))

#########################################################
# ACTUAL PROGRAM THAT USES GLOBAL VARIABLES AND FUNCTIONS
#########################################################

# Goes line by line in input file and computes individual contig coverage
with open(read_count, "rU") as rc_in:
    with open(output, "w") as out_f:
        for line in rc_in:
            line = line.strip()
            contig = line.split()[0]
            contig_lgth = lengths[contig]
            r_count = float(line.split()[1])
            coverage = cov(r_count, ave_read_lgth,contig_lgth)
            coverages.append(coverage)
            out_f.write(contig + "\t" + str(coverage) + "\n")

print "\n\033[1mMedian coverage\033[0m = ", median(coverages)
print "\n\033[1mAverage coverage\033[0m = ", ave_cov(coverages)
print "\n\033[1mRange\033[0m = ", range(coverages)[0], ",", range(coverages)[-1], "\n"
