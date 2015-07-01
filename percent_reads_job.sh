#!/bin/bash
#$ -N perc_read_map
#$ -M francois-olivier.gagnon-hebert.1@ulaval.ca
#$ -m beas
#$ -pe smp 8
#$ -l h_vmem=100G
#$ -l h_rt=10:00:00
#$ -cwd
#$ -S /bin/bash

time ./percent_reads_mapped.sh
