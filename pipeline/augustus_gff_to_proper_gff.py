#!/usr/bin/env python
from sys import argv, stderr

for line in open(argv[1]):
  if line[0] == '#': 
    print("#")
    continue

  atoms = line.split("\t")
  contig = atoms[0]
  prog = atoms[1]
  src = atoms[2]
  start = int(atoms[3])
  end = int(atoms[4])
  score = atoms[5]
  strand = atoms[6]
  frame = atoms[7]
  aug_attr = atoms[8][:-1]

  if src == 'gene':
    attr = "ID=" + aug_attr
  elif src == 'transcript':
    attr = "ID=%s;Parent=%s" % (aug_attr, aug_attr.split(".t")[0])
  elif src == 'CDS':
    aug_attr = aug_attr.split(";")[0].split()[1][1:-1]
    attr = "ID=%s.cds;Parent=%s" % (aug_attr, aug_attr)
  elif src == 'intron':
    aug_attr = aug_attr.split(";")[0].split()[1][1:-1]
    attr = "Parent=%s" % aug_attr
    

  print(f"{contig}\t{prog}\t{src}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr}")
  
  



 


  
