#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from re import compile
import numpy as np
from sys import argv, stderr
from collections import OrderedDict

f = SeqIO.parse(open(argv[1]), "fasta")

seq_d = {}
for rec in f:
  seq_d[rec.name] = str(rec.seq)

transcript_d = OrderedDict()

for line in open(argv[2]):
  if line[0] == '#': continue

  atoms = line.split("\t")
  contig = atoms[0]
  prog = atoms[1]
  src = atoms[2]
  start = int(atoms[3])
  end = int(atoms[4])
  strand = atoms[6]

  if src == 'CDS':
    transcript = atoms[-1].split()[1][1:-2]
    if strand == '+':
      transcript_d.setdefault(transcript, [])
      transcript_d[transcript].append((start, end, strand, seq_d[contig][start-1:end]))
    else:
      transcript_d.setdefault(transcript, [])
      transcript_d[transcript].append((start, end, strand, 
        str(Seq(seq_d[contig][start-1:end]).reverse_complement())))

for transcript in transcript_d:
  chunks = transcript_d[transcript]
  if chunks[0][2] == '+':
    chunks.sort()
  else:
    chunks.sort(reverse=True)
  
  seq = "".join([chunk[3] for chunk in chunks])
  print(">%s\n%s" % (transcript, seq))
    
  
  



 


  
