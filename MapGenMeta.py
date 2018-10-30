import os
import sys

##parse input output
##input meta
##aacc+"\t"+de+"\t"+source+"\t"+date+"\t"+strain+"\t"+ mol_type+"\t"+country+"\t"+ plasmid+ "\t" +isolate+"\n
out = open(sys.argv[3], 'w')
out.write("Acc\tFamily\tq.start\tq.end\tName\tOrganism\tDate\tStrain\tMolecularType\tCountry\tPlasmid\tIsolate\n")

meta_input = open(sys.argv[1], 'r')

g_input = open(sys.argv[2], 'r')
stat_out = open("stat_out", "w")

m = meta_input.readline()
acc = {}

while m:
    line = m.split("\t")
    line[2] = ' '.join(line[2].split()[:2])         ### refine source to species
    if (line[4] == ""):                             ##replace empty strain with isolate
        line[4] = line[8][:-2]
    acc[line[0]] = '\t'.join(line[1:])
    m = meta_input.readline()
meta_input.close()

#for key in rep_len.keys():
#    print "("+key+")"


query_id = ""
sub_id = ""
g = g_input.readline()
g = g_input.readline()
while g:
    g = g[:-1]
    l = g.split()
    l[0] = l[0][:-2]
    g = '\t'.join(l)
    if l[0] in acc:
        out.write(g+"\t"+acc[l[0]])
    else:
        out.write(g+"\t"+"NA\n")
    g =g_input.readline()
g_input.close()