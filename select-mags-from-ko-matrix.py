#!usr/bin/env python

import sys

fundict = {}
## Define the outfile prefix based on the name of the input file wich is the groupsing ids produced by ced and run through run-igraph-to-identify-groups.R
pre = sys.argv[1].split('-')[3]

fun_list = []
fun_names = []
#The sys.argv in the statement below should be the KEGG matrix "ALL-KEGG-O-table.txt"
for line in open(sys.argv[2], 'r'):
    x = line.strip().split('\t')
    if "MAG" == x[0]:
        fun_list.append(x)
    fundict[x[0]] = x[0:len(x)]
mag_group_list = []
mag_group_dict = {}

for line in open(sys.argv[1],'r'):
    x = line.strip().split('\t')
    if "V1" in x[0]:
        next
    else:
        mag_group_dict[x[0].strip('"')] = x[0:len(x)]
        mag_group_list.append(x[1])

unique_mag_group_list = set(sorted(mag_group_list))

# Collect only the groups that have between 2 and 20 members.
tight_mag_group_list = []
for i in unique_mag_group_list:
    if mag_group_list.count(i) >= 6 and mag_group_list.count(i) <= 50:
        print(i)
        tight_mag_group_list.append(i)

#Make a new dictionary of the MAGs that contain only the tight_mag_group_list
tight_mag_group_dict = {}
for i in mag_group_dict.keys():
#    print(i, mag_group_dict[i])
    if mag_group_dict[i][1] in tight_mag_group_list:
#        print(i, mag_group_dict[i])
        tight_mag_group_dict[i] = mag_group_dict[i]

print("collecting mags from the KO matrix for %d groups from %d total groups that contained between 10 and 50 members from %s network"%(len(tight_mag_group_list),len(unique_mag_group_list),pre))
print("also, just to be sure, we have %d MAGs in the tight groups which should match the number of lines in the KO matrix"%(len(tight_mag_group_dict.keys())))
### Writing a functional matrix for each of the groups in the network
for i in tight_mag_group_list:
    out = open('INDIVIDUAL-GROUP-FUNCTION-TABLES/'+pre+'-'+"group"+str(i)+'-'+'KO-matrix.txt', 'w')
    out.write('\t'.join(fun_list[0])+'\n')
    for key in tight_mag_group_dict.keys():
        if tight_mag_group_dict[key][1] == i:
            out.write('\t'.join(fundict[key])+'\n')

### Writing a complete functional matrix for all MAGs in the nutrient layer being examined KO.
outfile = open(pre+'-KO-matrix.txt','w')
outfile.write('\t'.join(fun_list[0])+'\n')
for i in tight_mag_group_dict.keys():
    outfile.write('\t'.join(fundict[i])+'\n')


### Writing a complete metadata matrix for all MAGs in the nutrient layer being examined.
meta_out = open(pre+'-KO-matrix-metadata.txt', 'w')
for i in open('mag-metadata-for-phyloseq.txt', 'r'):
    x = i.strip().split('\t')
    if "bin" in x[0]:
        meta_out.write('\t'.join(x)+'\t'+"net.group"+'\n')
    if x[0] in tight_mag_group_dict.keys():
#        print(mag_group_dict[x[0]][1])
        meta_out.write('\t'.join(x)+'\t'+str(tight_mag_group_dict[x[0]][1])+'\n')


        
### Writing a complete taxonomy matrix for all MAGs in the nutrient layer being examined
### This needs to be a two column list of each KO ID and the word "KEGG"
tax_out = open(pre+'-phyloseq-pseudo-tax.txt', 'w')
for i in fun_list[0]:
    tax_out.write(i+'\t'+"KEGG"+'\n')
