#!/usr/bin/env python
from sys import argv

if len(argv) > 2:
  min_splice_ratio = float(argv[2])
  min_introns = int(argv[3])
  max_intron_length = int(argv[4])
else:
  min_splice_ratio = 0.3
  min_introns = 10
  max_intron_length = 40

i = 0
for line in open(argv[1]).readlines():
  i += 1
  if line[0] != '#':
    atoms = line.split()

    start = int(atoms[1])# + 1
    end = int(atoms[2]) 

    intron_len = end - start + 1
    strand = atoms[3]
    #p = int(atoms[4])
    spliced_introns = int(atoms[5])
    fwd_strand_tot = int(atoms[6])
    rev_strand_tot = int(atoms[7])
    tot = int(atoms[8])

    if strand == '+' and fwd_strand_tot >= rev_strand_tot:
      sense_greater = True
    elif strand == '-' and rev_strand_tot >= fwd_strand_tot:
      sense_greater = True
    else:
      sense_greater = False

    if strand == '+':      
      splice_ratio = spliced_introns/(fwd_strand_tot + 0.01)
      #0.01 added to avoid division by 0 errors.
    elif strand == '-':
      splice_ratio = spliced_introns/(rev_strand_tot + 0.01)

    if sense_greater and splice_ratio >= min_splice_ratio and \
      spliced_introns >= min_introns and intron_len <= max_intron_length:

        print("%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\tID=intron%s" % (atoms[0],
                   "REALTRONS", "intron", start, end, splice_ratio, strand, ".", i))

# vim:sts=2:ts=2:sw=2:
