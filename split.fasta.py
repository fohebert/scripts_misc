#!/usr/local/bin/python

"""
\033[1mDESCRIPTION\033[0m
    Takes a fasta file and splits it into units of 'N'
    sequences.

\033[1mUSAGE\033[0m
    %program <input_fasta> <num_seq> <out_prefix>
"""

import sys
from Bio import SeqIO

try:
    in_fasta = sys.argv[1]
    num_seq = int(sys.argv[2])
    out_prefix = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

# - # - # - # - # - # -
# DEFINING FUNCTIONS #
# - # - # - # - # - # -

# Generator function that loads one record at a time and appends it 
# to a list. The process will be repeated until the list contains
# the desired number of sequences. This function is designed to
# deal with large fasta files (e.g. genomes or transcriptomes).
def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

# - # - # -
# PROGRAM #
# - # - # -

# Parsing the large fasta file provided as first argument
# in the command line ("in_fasta")
record_iter = SeqIO.parse(open(in_fasta),"fasta")

for i, batch in enumerate(batch_iterator(record_iter, num_seq)) :
    filename = "group_%i.fasta" % (i+1)
    handle = open(filename, "w")
    count = SeqIO.write(batch, handle, "fasta")
    handle.close()
    print "Wrote %i records to %s" % (count, filename)