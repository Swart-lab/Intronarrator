#!/usr/bin/env python
from Bio import SeqIO
from sys import argv, stderr

f = SeqIO.parse(open(argv[1]), "fasta")

for rec in f:
  print(">%s\n%s" % (rec.name, rec.seq.translate(table=4)))
