#!/usr/bin/env python
from Bio import SeqIO
from sys import argv, stderr

kind_d = {}
fh = open(argv[1])

while 1:
  has_intron = False

  l1 = fh.readline()
  if l1 == '': break
  l2 = fh.readline()
  l3 = fh.readline()
  
  while 1:
    aline = fh.readline()
    if aline == '\n':
      break

  if "intron" in l3:
    has_intron = True


  atoms1 = l1.split()
  atoms2 = l2.split()

  contig = atoms1[0].split(".")[0]
  tRNA_num = atoms1[0].split(".trna")[1]
  aa = atoms2[1]
  anticodon = atoms2[3]
  score = atoms2[-1]

  if anticodon != 'NNN':
    anticodon_loc = atoms2[6][1:-1].split("-")
    
    kind_d.setdefault((aa, anticodon), 0)
    kind_d[(aa, anticodon)] += 1

    prog = "tRNAscan-SE_2.0"
    inference=";inference=COORDINATES:profile:tRNAscan-SE:2.0"
    kind = "tRNA"
    a, b = int(atoms1[1].split("-")[0][1:]), int(atoms1[1].split("-")[1][:-1])
    if b > a:
      strand = '+'
      start, end = a, b
    elif b < a:
      strand = '-'
      start, end = b, a
    phase = "."
    gene_id = "ID=tRNA-%s(%s)_%s" % (aa, anticodon.lower(), kind_d[(aa, anticodon)])
    product = "product=tRNA-%s(%s)" % (aa, anticodon.lower())

    if int(anticodon_loc[0]) < int(anticodon_loc[1]):
      anticodon_str = ';anticodon="(pos:join(%s..%s),aa:%s,seq:%s)"' % (anticodon_loc[0], anticodon_loc[1], aa, anticodon)
    else:
      anticodon_str = ';anticodon="(pos:complement(%s..%s),aa:%s,seq:%s)"' % (anticodon_loc[0], anticodon_loc[1], aa, anticodon)

    note=''

    if has_intron:
      note=';note="possible intron: %s"' % l3.split()[2]

    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
      (contig, prog, kind, start, end, score, strand, phase, gene_id + ";" + product + inference + note))

  
