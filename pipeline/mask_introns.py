#!/usr/bin/env python
from Bio import SeqIO, Seq
from sys import argv

intron_d = {}

for line in open(argv[2]):
  atoms = line.split()
  intron_d.setdefault(atoms[0], []).append((int(atoms[3]), int(atoms[4])))

f = SeqIO.parse(open(argv[1]), "fasta")

for rec in f:
  seq_l = list(rec.seq)
  try:
    for intron in intron_d[rec.name]:
      for i in range(intron[0]-1, intron[1]):
        seq_l[i] = seq_l[i].lower()
  except KeyError:
    pass
  
  print(">%s\n%s" % (rec.name, ''.join(seq_l)))
