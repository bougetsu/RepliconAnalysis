import os
import sys

#os.system("sort -t $'\t' -k1,1r -k19,19n -k6,6n -k7,7n " + sys.argv[1] + " -o " + sys.argv[1]+"_filtered")
#a = "sort -t $'\t' -k1,1r -k19,19n -k6,6n -k7,7n " + sys.argv[1] + " -o " + sys.argv[1]+"_filtered"
#os.system(a)
#print a

stat_out = open("stat_out", "w")

hit_count = {}

query_id = ""
sub_id = ""
len = 0
start = 0
end = 0
family = ""
last_l = ""
in_filter = open(sys.argv[1]+"_filtered", 'r')
out_uniq = open(sys.argv[1]+"_uniq", 'w')
line = in_filter.readline()
out_uniq.write(line)
line = in_filter.readline()
while line:
    l = line.split("\t")
    if (query_id != l[0]):   ## new query
        if(last_l != ""):    ##new query, write last record
            out_uniq.write(last_l)
            hit_count[query_id] = hit_count[query_id] + 1
        last_l = line
        query_id = l[0]
        len = int(l[3])
        iden = float(l[2])
        start = int(l[5])
        end = int(l[6])
        family = l[18]
        print family+"\n"
        hit_count[query_id] = 0
    else:
        if (family == l[18]):
            if not (start > int(l[6]) or end < int(l[5])):   ##overlap
                if iden < float(l[2]):
                    len = int(l[3])
                    iden = float(l[2])
                    start = int(l[5])
                    end = int(l[6])
                    last_l = line
            else:      ##no overlap, same gene, new loc
                out_uniq.write(last_l)
                hit_count[query_id] = hit_count[query_id] + 1
                len = int(l[3])
                iden = float(l[2])
                start = int(l[5])
                end = int(l[6])
                family = l[18]
                last_l = line
        else: ## new gene
            out_uniq.write(last_l)
            hit_count[query_id] = hit_count[query_id] + 1
            len = int(l[3])
            iden = float(l[2])
            start = int(l[5])
            end = int(l[6])
            family = l[18]
            last_l = line
    line = in_filter.readline()
out_uniq.write(last_l)
in_filter.close()
out_uniq.close()


for key in hit_count.keys():
    stat_out.write(key+"\t"+str(hit_count[key])+"\n")

stat_out.close()