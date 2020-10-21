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

  clean_kind = kind.replace("_eukarya", "")
  clean_kind = clean_kind.replace("Protozoa_", "")
  gene_id = "%s-%s" % (clean_kind, kind_d[kind])

  if 'tRNA' in kind:
    continue

  ncRNA_kind = 'ncRNA'
  if "rRNA" in kind:
    ncRNA_kind = 'rRNA'
  elif kind[0] == "U":
    ncRNA_kind = 'snRNA'
  elif 'SRP' in kind:
    ncRNA_kind = 'SRP_RNA'
  elif 'MRP' in kind:
    ncRNA_kind = 'RNase_MRP_RNA'
  elif 'RNaseP' in kind:
    ncRNA_kind = 'RNase_P_RNA'
  

  print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s_gene" % 
    (contig, prog, "gene", start, end, score, strand, phase, gene_id))
  print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s" % 
    (contig, prog, ncRNA_kind, start, end, score, strand, phase, gene_id))

