#!/usr/bin/env python
from sys import argv

intron_d = {}
intron_d_starts = {}
intron_d_ends = {}

for line in open(argv[1]):
    atoms = line.split()
    s, e = int(atoms[3]), int(atoms[4])
    score = float(atoms[5])
    strand = atoms[6]
    id = atoms[-1]

    intron_d.setdefault(atoms[0], {})
    intron_d[atoms[0]][(s, e)] = (score, strand, id)
  
    #Populate dictionaries to detect overlapping alternatively spliced introns
    #intron_d_starts[(atoms[0], s)].append((e, score))
    for i in range(s, e +1):
      intron_d_starts.setdefault((atoms[0], i), [])
      intron_d_starts[(atoms[0], i)].append((e, score))

    #intron_d_ends[(atoms[0], e)].append((s, score))
    for i in range(s, e +1):
      intron_d_ends.setdefault((atoms[0], i), [])
      intron_d_ends[(atoms[0], i)].append((s, score))

new_d = {}
for contig in intron_d:
  for intron_s, intron_e in intron_d[contig]:
    intron_d_e_len = len(intron_d_ends[(contig, intron_e)])
    intron_d_s_len = len(intron_d_starts[(contig, intron_s)])

    new_d.setdefault(contig, {})
    if intron_d_s_len == 1 and intron_d_e_len == 1:
      new_d[contig][(intron_s, intron_e)] = intron_d[contig][(intron_s, intron_e)]
    else:
      if intron_d_e_len > 1:
        e_max_score = max([score for start, score in intron_d_ends[(contig, intron_e)]])
        #print("e_max_score", e_max_score, intron_d[contig][intron_s, intron_e][0])

        #record only the most efficiently spliced alternatively excised intron
        if intron_d[contig][intron_s, intron_e][0] == e_max_score:
          new_d.setdefault(contig, {})
          new_d[contig][(intron_s, intron_e)] = intron_d[contig][(intron_s, intron_e)]

      elif intron_d_s_len > 1:
        s_max_score = max([score for start, score in intron_d_starts[(contig, intron_s)]])
        if intron_d[contig][intron_s, intron_e][0] == s_max_score:
          new_d[contig][(intron_s, intron_e)] = intron_d[contig][(intron_s, intron_e)]
       
for contig in new_d:
  if new_d[contig] != {}:
    for intron_s, intron_e in new_d[contig]:
      score, strand, id = new_d[contig][(intron_s, intron_e)]
      print("%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s" % (contig, "REALTRONS", "intron", intron_s, intron_e, score, strand, '.', id))
