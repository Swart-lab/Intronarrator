#!/usr/bin/env python
from Bio import SeqIO
from sys import argv, stderr

kind_d = {}
for line in open(argv[1]).readlines():  
  atoms = line.split()
  contig = atoms[0]
  kind = atoms[2]

  kind_d.setdefault(atoms[2], 0)
  kind_d[atoms[2]] += 1

  prog = "Infernal"
  score = atoms[-3]

  a, b = int(atoms[7]), int(atoms[8])
  if b > a:
    strand = '+'
    start, end = a, b
  elif b < a:
    strand = '-'
    start, end = b, a
  phase = "."
  gene_id = "%s-%s" % (kind, kind_d[kind])

  print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
    (contig, prog, 'ncRNA', start, end, score, strand, phase, gene_id))

