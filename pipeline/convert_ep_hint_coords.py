#!/usr/bin/env python
from Bio import SeqIO, Seq
from re import compile
import numpy as np
from sys import argv, stderr

p = compile("[actg]")
f = SeqIO.parse(open(argv[1]), "fasta")
#lowercase masked intron contigs

ep_hint_length = 10

contig_coords = {}
for rec in f:
  aseq = str(rec.seq)
  strip_seq = p.sub('', aseq)
  coords = np.zeros((len(aseq) + 1,), dtype=int)

  newseq = ""
  coord_d = {}
  i, j = 0, 0

  for b in aseq:
    i += 1
    if b.islower():
      j += 1
    else:
      coords[i] = i - j
    
  contig_coords[rec.name] = coords

for line in open(argv[2]):
  atoms = line.split()
  contig = atoms[0]
  prog = atoms[1]
  src = atoms[2]
  
  try:
    s, e = contig_coords[contig][int(atoms[3])], contig_coords[contig][int(atoms[4])]
    score = atoms[5]
    strand = atoms[6]
    note = atoms[-1]

    if e - s + 1 == ep_hint_length:
      print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, src, s, e, score, strand, ".", note))

  except IndexError:
    stderr.write("EP hint out of bounds: %s %s %s" % (contig, atoms[3], atoms[4]))
