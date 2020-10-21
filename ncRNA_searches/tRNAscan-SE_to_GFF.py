#!/usr/bin/env python
from Bio import SeqIO
from sys import argv, stderr

kind_d = {}
for line in open(argv[1]).readlines()[3:]:  
  atoms = line.split()
  contig = atoms[0]
  tRNA_num = int(atoms[1])
  aa = atoms[4]
  anticodon = atoms[5]
  
  kind_d.setdefault((aa, anticodon), 0)
  kind_d[(aa, anticodon)] += 1

  prog = "tRNAscan-SE_2.0"
  kind = "tRNA"
  if atoms[-1] != 'pseudo':
    score = float(atoms[-1])
  else:
    score = 0
  a, b = int(atoms[2]), int(atoms[3])
  if b > a:
    strand = '+'
    start, end = a, b
  elif b < a:
    strand = '-'
    start, end = b, a
  phase = "."
  gene_id = "ID=tRNA-%s(%s)-%s" % (aa, anticodon, kind_d[(aa, anticodon)])

  print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
    (contig, prog, kind, start, end, score, strand, phase, gene_id))

  if int(atoms[6]) > 0:
    if b > a: 
      print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
          (contig, prog, "intron", atoms[6], atoms[7], score, strand, phase, gene_id))
    elif b < a:
      print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
        (contig, prog, "intron", atoms[7], atoms[6], score, strand, phase, gene_id))
  
