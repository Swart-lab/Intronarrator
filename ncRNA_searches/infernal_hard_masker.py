#!/usr/bin/env python
from Bio import SeqIO
from sys import argv, stderr

f = SeqIO.parse(open(argv[2]), "fasta")

ncRNA_d = {}
for line in open(argv[1]):  
  if line[0] != '#':
    atoms = line.split()
    e_val = float(atoms[15])
    if e_val < 1e-6:
      ncRNA_d.setdefault(atoms[0], []).append([int(atoms[7]), int(atoms[8])])
      #print(line[:-1])

for rec in f:
  seq_l = list(rec.seq)
  try:
    for ncRNA in ncRNA_d[rec.name]:
      for i in range(ncRNA[0]-1, ncRNA[1]):
        seq_l[i] = 'N'
  except KeyError:
    pass
  
  print(">%s\n%s" % (rec.name, ''.join(seq_l)))
