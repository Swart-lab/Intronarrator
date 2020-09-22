#!/usr/bin/env python
from Bio import SeqIO, Seq
from re import compile
import numpy as np
from sys import argv

p = compile("[actg]")
f = SeqIO.parse(open(argv[1]), "fasta")

for rec in f:
  aseq = str(rec.seq)
  strip_seq = p.sub('', aseq)
  coords = np.zeros((len(strip_seq) + 1,), dtype=int)

  newseq = []
  coord_d = {}
  i, j = 0, 0

  for b in aseq:
    i += 1
    if b.islower():
      j += 1
    else:
      newseq.append(b)
      #coords[i - j] = i
      #coord_d[i - j] = i, b
    
  print(">%s\n%s" % (rec.name, "".join(newseq)))
  #print(coords)
  #print(len(newseq), len(rec.seq))

