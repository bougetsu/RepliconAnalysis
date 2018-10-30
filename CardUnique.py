import os
import sys


stat_out = open("Card_overlap_difgene", "w")

hit_count = {}

query_id = ""
sub_id = ""
len = 0
qstart = 0
qend = 0
sstart = 0
send = 0
last_l = ""
iden = 0.0

slen = 0

reverse = False
in_filter = open(sys.argv[1]+"_sorted", 'r')
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
        print "new query" + l[0]+ "\n"
        last_l = line
        query_id = l[0]
        sub_id = l[1]
        len = int(l[3])
        iden = float(l[2])
        if( int(l[5]) > int(l[6]) ):
            reverse = True
        else:
            reverse = False
        qstart = int(l[5])
        qend = int(l[6])
        sstart = int(l[7])
        send = int(l[8])
        hit_count[query_id] = 0
        print l[0]+ "\t" + str(qstart) + "\t" + str(qend) + "\t" +  "Reverse"+ str(reverse)+  "\n"
    else:
        print "Same query new line " + l[5] + "\t" + l[6] + "\n"
        if ( (reverse == False) and (qstart > int(l[6]) or qend < int(l[5])) ) or ( (reverse == True) and (qstart < int(l[6]) or qend > int(l[5])) ):  #####no overlap
            print l[0] + "NO Overlap \n"
            out_uniq.write(last_l)
            hit_count[query_id] = hit_count[query_id] + 1
            sub_id = l[1]
            len = int(l[3])
            iden = float(l[2])
            qstart = int(l[5])
            qend = int(l[6])
            last_l = line
        else:
            ####overlap length
            slen = min(len , int(l[3]))
            if (reverse == False):
                overlen = min( (qend - int(l[5])), (int(l[6]) - qstart) )
            else:
                overlen = min(( -qend + int(l[5])), ( -int(l[6]) + qstart))

            print l[0] + " Overlap Overlen " + str(overlen) +  " slen " + str(slen) + " \n"
            ###small or larger overlap
            if ( float(overlen)/ float(slen) > 0.5 ): ####larger overlap at same loc, retain the highest identity
                print "Large overlap \n"
                if (iden < float(l[2])):
                    last_l = line
            else:
                print "small overlap \n"
                if (sub_id == l[1]):
                    if (reverse == False):
                        l[5] = str(min(int(l[5]), qstart))
                        l[6] = str(max(int(l[6]), qend))
                        l[7] = str(min(int(l[7]), sstart))
                        l[8] = str(max(int(l[8]), send))
                    else:
                        l[5] = str(max(int(l[5]), qstart))
                        l[6] = str(min(int(l[6]), qend))
                        l[7] = str(max(int(l[7]), sstart))
                        l[8] = str(min(int(l[8]), send))
                    last_l = "\t".join(l)
                else:
                    stat_out.write(last_l)
                    stat_out.write(line)
                    stat_out.write("##############\n")
                    out_uniq.write(last_l)
                    last_l = line
                    query_id = l[0]
                    sub_id = l[1]
                    len = int(l[3])
                    iden = float(l[2])
                    iden = float(l[2])
                    if (int(l[5]) > int(l[6])):
                        reverse = True
                    else:
                        reverse = False
                    qstart = int(l[5])
                    qend = int(l[6])
                    sstart = int(l[7])
                    qend = int(l[8])
    line = in_filter.readline()
out_uniq.write(last_l)
in_filter.close()
out_uniq.close()

stat_out.close()