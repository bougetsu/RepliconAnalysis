import os
import sys


out = open(sys.argv[1]+"ResName", 'w')
out.write("Acc\tLabel\tFeature\tq. start\tq. end\tName\tOrganism\n")

rep_input = open(sys.argv[1], 'r')

line = rep_input.readline()

while line:
    l = line.split("\t")
    print "\t".join(l)
    if(l[1] == "Res"):
        res = l[2].split("|")
        l[2] = res[-1]
    out.write("\t".join(l))
    line = rep_input.readline()
rep_input.close()
