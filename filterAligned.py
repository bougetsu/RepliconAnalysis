import os
import sys

##parse input output
##input rep_info.tab
##Gene/ORF	Species	Plasmid	GenBank_no.	Position	Size	Rep-family	Protein	Domain I	Domain II	seq_id
out = open(sys.argv[3], 'w')
out.write("quer_ id\tsubject_id\t% identity\talignment length\t% query coverage per subject\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tseq_id\tGene/ORF\tSpecies\tPlasmid\tGenBank_no.\tPosition\tSize\tRep-family\tProtein\tDomain_I\tDomain_II\n")

rep_input = open(sys.argv[1], 'r')

blast_input = open(sys.argv[2], 'r')
stat_out = open("stat_out", "w")

rep = rep_input.readline()
rep_info = {}
rep_len = {}
hit_count = {}

while rep:
    rep = rep[:-2]
    info = rep.split("\t")
    rep_info[info[0]] = rep
    rep_len[info[0]] =  info[6]
    rep = rep_input.readline()
rep_input.close()

#for key in rep_len.keys():
#    print "("+key+")"


query_id = ""
sub_id = ""
blast = blast_input.readline()
while blast:
    blast = blast[:-2]
    if not (blast.startswith("#")):
        line = blast.split('\t')
        if (query_id != line[0]):
            query_id = line[0]
            hit_count[query_id] = 0
        sub_id = line[1]
        align_len = line[3]
        if (float(line[2]) >= 80) and ( float(align_len)/ float(rep_len[sub_id]) >= 0.6):
            out.write(blast+"\t"+rep_info[sub_id]+"\n")

    blast = blast_input.readline()
blast_input.close()

