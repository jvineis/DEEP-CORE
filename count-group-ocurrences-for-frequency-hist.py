#!/usr/bin/env python
import sys
from collections import Counter

outfile = open(sys.argv[2], 'w')
outfile.write("group"+'\t'+"count"+'\t'+"layer"+'\n')
depth_group = []
for i in open(sys.argv[1], 'r'):
    x = i.strip().split('\t')
    depth_group.append(x[2])

depth_group_count = Counter(depth_group)

for key in depth_group_count.keys():
    if depth_group_count[key] == 1:
        next
    else:
        outfile.write(str(key)+'\t'+str(depth_group_count[key])+'\t'+key.split("_")[0]+'\n')
    
