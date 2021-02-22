#!/usr/bin/env python

import sys

outfile = open(sys.argv[2], 'w')
outfile.write("name" +'\t'+"contigs_db_path" + '\n')

for i in open(sys.argv[1], 'r'):
    x = i.strip()
    print(x)
    outfile.write(x+'\t'+x+".db"+'\n')

