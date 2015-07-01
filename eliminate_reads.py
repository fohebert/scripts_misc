#!/usr/bin/python

import sys
import re
from Bio import SeqIO

# Takes a .qual file and parse it to eliminate certain reads that
# are unwanted (e.g. reads containing transposon sequences) based
# on a text file containing the unwanted names.

input_qual = sys.argv[1]
input_fasta = sys.argv[2]
unwanted_reads = sys.argv[3] # One unwanted read per line
output_file = sys.argv[4]    

# Returns a sfastq file without the unwanted names

unwanted_names = set()

with open(unwanted_reads, "r") as f:
	for line in f:
		line = line.strip()
		unwanted_names.add(line)

reads = SeqIO.to_dict(SeqIO.parse(open(input_fasta), "fasta"))

with open(output_file, "w") as f:
	for rec in SeqIO.parse(open(input_qual), "qual"):
		reads[rec.id].letter_annotations["phred_quality"] = \
		rec.letter_annotations["phred_quality"]		
		if rec.id not in unwanted_names:
			f.write(reads[rec.id].format("fastq"))
