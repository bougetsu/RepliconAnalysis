import os
import sys

##parse input output
##input rep_info.tab
##Gene/ORF	Species	Plasmid	GenBank_no.	Position	Size	Rep-family	Protein	Domain I	Domain II	seq_id
out = open("AllIdGenome", 'w')

input = open(sys.argv[1], 'r')

line = input.readline()

while line:
    if (line.startswith("# Query")):
        l = line.split(' ')
        id = l[2]
        id = id[:-2]
        out.write(id+"\n")
    line = input.readline()
input.close()
out.close()