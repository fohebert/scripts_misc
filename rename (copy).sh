#!/usr/bin/env bash

for i in `ls -1 *.sfastq | sed 's/lb_[a-z]*\.pl_illumina\.sm_//g' | sed 's/\.sfastq//g' | sort -u |`; do ./rename_reads_sfastq_illumina.py $i_merged.sfastq $i $i_new_name.sfastq; done
