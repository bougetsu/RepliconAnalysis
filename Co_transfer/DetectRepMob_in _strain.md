###Co-transfer of plasmids

####Steps

1. Data: genome downded from ncbi 2018.03
2. Process of Data: __remove redudant records!!!__
3. Using Mob_suit to get MOB, 


####Remove redudant records

1. From meta data, find same or similar strains

```r
acas <- read.table("ncbi-genomes-2018-03-19/acc_fasta2", header = T, sep = "\t", stringsAsFactors = F)

dup_Strain <- acas %>% left_join(alist, by = "Assembly") %>%
  select(Assembly, Organism, Strain) %>%
  group_by(Organism, Strain) %>%
  distinct(Assembly) %>%
  summarise(n = n()) %>%
  filter(n != 1) %>%
  select(Organism, Strain)
dim(dup_Strain)
#58 2
dup_assembly <- dup_Strain %>% left_join(alist, by = c("Organism", "Strain")) %>%
  select(Assembly, Organism, Strain) %>%
  group_by(Organism, Strain) %>%
  filter(row_number()!=1)

write.table(dup_assembly, file = "Co_transfer/dup_assembly.tsv", sep = "\t", row.names = F, quote = F)
```

2. After remove those redundant strains, look into details about similar strain/isolates

```r
acas %>% left_join(alist, by = "Assembly") %>%
  select(Assembly, Organism, Strain) %>%
  group_by(Organism, Strain) %>%
  distinct(Assembly) %>%
  filter(! Assembly %in% dup_assembly$Assembly) %>%
  arrange(Organism, Strain) %>%
  write.table(file = "Co_transfer/check_assembly_list.tsv", sep = "\t", row.names = F, quote = F)
```
