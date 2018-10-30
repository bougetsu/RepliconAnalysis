import os
import sys

##parse input output
##input rep_info.tab
##Gene/ORF	Species	Plasmid	GenBank_no.	Position	Size	Rep-family	Protein	Domain I	Domain II	seq_id
out = open(sys.argv[2], 'w')
out.write("quer_ id\tsubject_id\t% identity\talignment length\t% query coverage per subject\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tseq_id\tGene/ORF\tSpecies\tPlasmid\tGenBank_no.\tPosition\tSize\tRep-family\tProtein\tDomain_I\tDomain_II\n")

blast_input = open(sys.argv[1], 'r')


query_id = ""
sub_id = ""
blast = blast_input.readline()
while blast:
    blast = blast[:-2]
    if not (blast.startswith("#")):
        line = blast.split('\t')
        if (query_id != line[0]):
            query_id = line[0]
        sub_id = line[1]
        align_len = line[3]
        if (float(line[2]) >= 80) :
            out.write(blast+"\n")
    blast = blast_input.readline()
blast_input.close()

