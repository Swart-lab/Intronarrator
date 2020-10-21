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
  extra_id = ''

  if "rRNA" in kind:
    ncRNA_kind = 'rRNA'
  elif kind[0] == "U":
    extra_id = ';ncRNA_class="snRNA"'
  elif 'SRP' in kind:
    extra_id =  ';ncRNA_class="SRP_RNA"'
  elif 'MRP' in kind:
    extra_id = ';ncRNA_class="RNase_MRP_RNA"'
  elif 'RNaseP' in kind:
    extra_id = ';ncRNA_class="RNase_P_RNA"'
  

  print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s_gene%s" % 
    (contig, prog, "gene", start, end, score, strand, phase, gene_id, extra_id))
  print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s%s" % 
    (contig, prog, ncRNA_kind, start, end, score, strand, phase, gene_id, extra_id))

