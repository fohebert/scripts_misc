#!/bin/bash
#$ -N read-lgth
#$ -M francois-olivier.gagnon-hebert.1@ulaval.ca
#$ -m eas
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=10:00:00
#$ -cwd
#$ -S /bin/bash

for read_file in `ls -1 *.R1* | sed 's/\.R1\.trim\.paired\.fastq//g'`; do
    cat ${read_file}.R1.trim.paired.fastq ${read_file}.R2.trim.paired.fastq | \
    awk 'NR%4 == 2 {s += length($1); t++} END {print s/t}' >${read_file}.ave.read-lgth;
done

cat *.ave.read-lgth >all-samples.read-length.txt

rm *.ave.read-lgth
