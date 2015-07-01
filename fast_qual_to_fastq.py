#!/usr/bin/python

from Bio import SeqIO
import sys
import re

# Converts a FASTA file and a QUAL file from Roche 454 sequencing in a FASTQ file.

quality_file = sys.argv[1] # 454 output quality file <filename>.qual
fasta_file = sys.argv[2] # 454 output fasta file : <filename>.fna
output_file = sys.argv[3] # Out-put file with FASTQ format

reads = SeqIO.to_dict(SeqIO.parse(open(fasta_file), "fasta"))

with open(output_file, "w") as f:
    for rec in SeqIO.parse(open(quality_file), "qual"):
        reads[rec.id].letter_annotations["phred_quality"] = rec.letter_annotations["phred_quality"]
        f.write(reads[rec.id].format("fastq"))
