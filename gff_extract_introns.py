#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from re import compile
import numpy as np
from sys import argv, stderr

f = SeqIO.parse(open(argv[1]), "fasta")

seq_d = {}
for rec in f:
  seq_d[rec.name] = str(rec.seq)

for line in open(argv[2]):
  if line[0] == '#': continue

  atoms = line.split()
  contig = atoms[0]
  prog = atoms[1]
  src = atoms[2]
  start = int(atoms[3])
  end = int(atoms[4])
  strand = atoms[6]

  if src == 'intron':
    if strand == '+':
      intron_seq = seq_d[contig][start-1:end]
      if intron_seq[:2].lower() != 'gt':
      #if intron_seq[-2:].lower() != 'ag':
        print(intron_seq, contig, atoms[-1], strand, start, end)
    else:
      intron_seq = Seq(seq_d[contig][start-1:end]).reverse_complement()
      if intron_seq[:2].lower() != 'gt':
      #if intron_seq[-2:].lower() != 'ag':
        print(intron_seq, contig, atoms[-1], strand, start, end)




  
