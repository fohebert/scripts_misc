#!/bin/bash
# Read alignment for quality assessment - RNA-seq metrics - % reads mapped to assembly

/prg/trinityrnaseq/trinityrnaseq_r20140717/util/align_and_estimate_abundance.pl --transcripts ../05_trinity_output/Trinity.fasta --seqType fq --left ./FOGrand_S1_L001_R1.paired.fastq --right FOGrand_S1_L001_R2.paired.fastq --est_metho RSEM --aln_method bowtie2 --trinity_mode --prep_reference
