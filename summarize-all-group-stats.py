#!/usr/bin/env python

import sys

outfile = open(sys.argv[2], 'w')
outfile.write("depth"+'\t'+"group_size"+'\t'+"number_of_groups"+'\t'+"different_genera"+'\t'+"different_family"+'\t'+"different_core"+'\t'+"different_depth"+'\t'+"different_creek"+'\n')
sizes = []
table_dict = {}
def summarize_depth_and_group(depth, size):
    genus_result = []
    family_result = []
    core_result = []
    depth_result = []
    creek_result = []

    results = []

    for key in table_dict.keys():
        if table_dict[key][3] == size and table_dict[key][0] == depth:
            print(table_dict[key])
            if int(table_dict[key][2]) == 1:
                genus_result.append(0)
            if int(table_dict[key][2]) > 1:
                genus_result.append(1)
            if int(table_dict[key][4]) == 1:
                family_result.append(0)
            if int(table_dict[key][4]) > 1:
                family_result.append(1)
            if int(table_dict[key][5]) == 1:
                core_result.append(0)
            if int(table_dict[key][5]) > 1:
                core_result.append(1)
            if int(table_dict[key][6]) == 1:
                depth_result.append(0)
            if int(table_dict[key][6]) > 1:
                depth_result.append(1)
            if int(table_dict[key][7]) == 1:
                creek_result.append(0)
            if int(table_dict[key][7]) > 1:
                creek_result.append(1)

    p_genus = sum(genus_result)/len(genus_result)
    p_family = sum(family_result)/len(family_result)
    p_core = sum(core_result)/len(core_result)
    p_depth = sum(depth_result)/len(depth_result)
    p_creek = sum(creek_result)/len(creek_result)
    results.append([depth,size,len(genus_result),p_genus, p_family, p_core, p_depth, p_creek])
    print(results)
    return(results)

            
for line in open(sys.argv[1], 'r'):
    x = line.strip().split('\t')
    if x[0] == "depth":
        next
    else:
        x_key = x[0]+'-'+x[1]
        table_dict[x_key] = x[0:len(x)]
        sizes.append(x[3])

u_sizes = set(sizes)

print(u_sizes)

for group_size in u_sizes:
    print(group_size)
    sud = summarize_depth_and_group("deep", group_size)
    outfile.write('\t'.join(str(x) for x in sud[0])+'\n')
for group_size in u_sizes:
    summ = summarize_depth_and_group("mid", group_size)
    outfile.write('\t'.join(str(x) for x in summ[0])+'\n')
for group_size in u_sizes:
    sus = summarize_depth_and_group("shallow", group_size)
    outfile.write('\t'.join(str(x) for x in sus[0])+'\n')
