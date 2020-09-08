#!/usr/bin/env python
from Bio import SeqIO
from sys import argv, stderr

#f = SeqIO.parse(open(argv[1]), "fasta")

#for rec in f:
#  print(">%s\n%s" % (rec.name, rec.seq.translate(table=4)))

not_allowed = ['bacteria', 'archaea', 'microsporidia', 'Plant', 'ACEA',
'etazoa', 'ictyostelium', 'RMST_9']

for line in open(argv[1]):  
  if line[0] != '#':
    atoms = line.split()
    e_val = float(atoms[15])
    skip = False
    if e_val < 1e-6:
      for nn in not_allowed:
        if nn in atoms[2]:
          skip = True

      if not skip: 
        print(line[:-1])

