#!/usr/bin/python

"""
DESCRIPTION
        Takes a fasta file and splits it into files of 
        <numer_entries_per_file> (integer) entries per file.

USAGE
        %program <input_FASTA_file> <number_entries_per_file>
        
        Note : the output files are named "scaffold_#.fasta"
               it is possible to change it from the source
               code.

"""

import sys
from Bio import SeqIO

try:
    in_file = sys.argv[1]
    number_entries = int(sys.argv[2])
except:
    print __doc__
    sys.exit(0)

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

record_iter = SeqIO.parse(open(in_file), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, number_entries)):
    filename = "scaffold_%i.fasta" % (i+1)
    handle = open(filename, "w")
    count = SeqIO.write(batch, handle, "fasta")
print "Done creating files"
