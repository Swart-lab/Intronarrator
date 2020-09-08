#!/usr/bin/env python
from Bio import SeqIO, Seq
from re import compile
import numpy as np
from sys import argv
import sys

p = compile("[actg]")
f = SeqIO.parse(open(argv[1]), "fasta")

contig_coords = {}
for rec in f:
  aseq = str(rec.seq)
  strip_seq = p.sub('', aseq)
  coords = np.zeros((len(strip_seq) + 1,), dtype=int)

  newseq = ""
  coord_d = {}
  i, j = 0, 0

  for b in aseq:
    i += 1
    if b.islower():
      j += 1
    else:
      #newseq += b
      #coord_d[i - j] = i, b
      coords[i -j] = i
    
  contig_coords[rec.name] = coords

realtron_d = {}
for line in open(argv[3]):
  atoms = line.split()
  realtron_d.setdefault(atoms[0], []).append((int(atoms[3]), int(atoms[4]), atoms[5], atoms[6]))
  
for line in open(argv[2]):
  if line[0] != '#':
    atoms = line.split("\t")
    contig = atoms[0]

    #if contig == "contig_83":
    try:
      prog = atoms[1]
    except KeyError:
      sys.stderr.write("curr line %s\n" % line)
    
    src = atoms[2]

    #try:
    s, e = contig_coords[contig][int(atoms[3])], contig_coords[contig][int(atoms[4])]
    score = atoms[5]
    strand = atoms[6]
    note = atoms[-1][:-1]

    wrong_strand = False

    #if r'"g6"' in note or r'"g70"' in note:
    #Debugging code above

    if src == 'CDS':
      #if "contig_83.g10.t1" in note: print("here")
      cds_coords = [s, e]
      try:
        for intron in realtron_d[contig]:
          if intron[0] >= s and intron[1] <= e:
            if intron[3] == strand:
              #check to make sure realtrons intron orientation corresponds to
              #AUGUSTUS prediction, otherwise remove bad AUGUSTUS gene prediction
              #if "contig_83.g10.t1" in note: print("here2")

              cds_coords.append(intron[0] - 1 )
              cds_coords.append(intron[1] + 1)
            else:
              wrong_strand = True
      except KeyError:
        #ignore cases where no intron is present in RNA-seq data for a contig/scaffold
        pass

      cds_coords.sort()         

      if not wrong_strand:

        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, "gene", s, e,
  score, strand, ".", note.split()[-1][1:-1]))
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, "transcript", s, e,
  score, strand, ".", note.split()[1][1:-2]))
        
        if strand == '+':
          phase = 0
          c_s, c_e = cds_coords[0], cds_coords[1]
          print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, src, c_s, c_e, score, strand, 0, note))

          if len(cds_coords) > 2:
            for i in range(2, len(cds_coords), 2):
              c_s, c_e = cds_coords[i], cds_coords[i+1]
              rem = ((3-phase) + cds_coords[i-1] - cds_coords[i-2] + 1) %3
              
              if rem == 0: 
                phase = 0
              elif rem == 1:
                phase = 2
              else:
                phase = 1

              print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, src, c_s,
                c_e, score, strand, phase, note))

        else:
          phase = 0
          c_s, c_e = cds_coords[-2], cds_coords[-1]
          print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, src, c_s, c_e, score, strand, 0, note))

          if len(cds_coords) > 2:
            for i in range(len(cds_coords) -3, 0, -2):
              #sys.stderr.write("i: %s\n %s" % (i, cds_coords))
              #Debugging code above - antisense is hard!

              c_s, c_e = cds_coords[i-1], cds_coords[i]
              seg_length = abs(cds_coords[i+1] - cds_coords[i+2]) + 1
              rem = ((3-phase) + seg_length) %3
              
              if rem == 0: 
                phase = 0
              elif rem == 1:
                phase = 2
              else:
                phase = 1

              print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, prog, src, c_s,
                c_e, score, strand, phase, note))

        if len(cds_coords) > 2:
          for i in range(1, len(cds_coords) -1, 2):
            i_s, i_e = cds_coords[i] + 1 , cds_coords[i+1] -1
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, "INTRONARRATOR", "intron", i_s, i_e, score, strand, ".", note))
            #TODO: get score from INTRONARRATOR

        print("#") 

  #except IndexError:
  #  sys.stderr.write(line)
