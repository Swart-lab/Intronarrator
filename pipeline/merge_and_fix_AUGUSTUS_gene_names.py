#!/usr/bin/env python

#USAGE: augustus_merge_and_fix_names.py "*.gff" > merged.gff
#NOTE the quotes, which are necessary to prevent shell expansion of the "*"

import sys
from sys import argv
from glob import glob

fns = [fn for fn in glob(argv[1])]

max_fn = max([int(fn.split(".")[-2]) for fn in fns])

prefix = ".".join(fns[0].split(".")[:-2])
suffix = fns[0].split(".")[-1]

for i in range(1, max_fn+1):
  fh = open(prefix + "." + str(i) + "." + suffix)

  for line in fh:
    #if 'start gene' in line or 'end gene' in line:
    #  print(line[:-1])

    if line[0] != '#':

      atoms = line.split("\t")

      try:
        if atoms[2] == 'gene':
          print("\t".join(atoms[:-1]) + "\t" + atoms[0] + "." + atoms[-1][:-1])
        elif atoms[2] == 'transcript':
          print("\t".join(atoms[:-1]) + "\t" + atoms[0] + "." + atoms[-1][:-1])
        elif atoms[2] == 'CDS':
          transcript = atoms[0] + "." + atoms[-1].split(";")[0].split()[1][1:-1]
          gene = atoms[0] + "." + transcript.split(".")[1]
          print("\t".join(atoms[:-1]) + "\t" + 'transcript_id "%s"; gene_id "%s"' % (transcript, gene))
      except IndexError:
        sys.stderr.write("line: " + line)

      

# vim:sts=2:ts=2:sw=2
