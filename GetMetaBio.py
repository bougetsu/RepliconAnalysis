from Bio import Entrez
from Bio import SeqIO
Entrez.email = "zhcong.pku@gmail.com"


with open('AllIdGenome') as f:
    id_list = f.read().splitlines()


id_list
id=",".join(id_list)

fetch_handle = Entrez.efetch(db="nucleotide", retmode="xml", id=id)
data = Entrez.read(fetch_handle)
fetch_handle.close()

out_seq = open("All_geno_meta.tab", "w")


for record in data:
	de = record['GBSeq_definition']
	source = record['GBSeq_source']
	acc = record['GBSeq_primary-accession']
	date = record['GBSeq_create-date']
	strain = ""
	isolate = ""
	country = ""
	plasmid = ""
	mol_type = ""
	orgn = ""
	for i in record['GBSeq_feature-table'][0]['GBFeature_quals']:
		if (i.values()[0] == "strain"):
			strain = i.values()[1]
		elif (i.values()[0] == "isolate"):
			isolate = i.values()[1]
		elif (i.values()[0] == "country"):
			country = i.values()[1]
		elif (i.values()[0] == "plasmid"):
			plasmid = i.values()[1]
		elif (i.values()[0] == "mol_type"):
			mol_type = i.values()[1]
		elif (i.values()[0] == "organism"):
			orgn = i.values()[1]
	print acc, de, orgn, source, mol_type, orgn,date
	out_seq.write(acc+"\t"+de+"\t"+source+"\t"+date+"\t"+strain+"\t"+ mol_type+"\t" +isolate+"\t"+country+"\t"+ plasmid+"\n")
out_seq.close()