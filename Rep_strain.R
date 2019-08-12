setwd("~/Work/project/data/")



###assembly
rep <- read.table("Assembly_Rep_nodup.xls", sep = "\t", header = T, stringsAsFactors = F)
dat <- read.table("enter_geno_family_meta.tab", header = T, stringsAsFactors = F, sep = "\t",row.names = NULL, quote = "")
alist <- read.csv("TableEnterococci.csv", header = T, stringsAsFactors = F)
alist$X.Organism.Name <- gsub("(\\w*) (\\w*) (.*)", "\\1 \\2",alist$X.Organism.Name)
alist$X.Organism.Name <- gsub("(\\w*) (\\w*)\\. (.*)", "\\1 \\2\\.",alist$X.Organism.Name)
colnames(alist)[1] <- "Organism"

dupgca <- read.table("Dup_GCA")
colnames(dupgca) <- c("Assembly")
acas <- read.table("ncbi-genomes-2018-03-19/acc_fasta2", header = T, sep = "\t", stringsAsFactors = F)
colnames(acas)[2] <- "Acc"
acas$Acc <- gsub("\\.[1-9]","",acas$Acc)


##dup accession

dupacc <- merge(dupgca, acas, by = "Assembly")

dupgca %>% left_join(acas, by = "Assembly")

##remove dupped acc and ass in rep blast data
##remove duplicated CGA and corresponding ACC FROM DATA
idx_gca <- which( !(alist$Assembly %in% dupgca[,1]) )

alist_nodup <- alist[idx_gca, ]

##find new dup
alist_nodup[which(alist_nodup$Strain  %in% alist_nodup[which(duplicated(alist_nodup$Strain)),3]), c(1, 3, 6)]

##merge with rep results

acas <- acas[which(acas$Assembly %in% alist_nodup$Assembly),]

dataa <- left_join(dat, acas, by = "Acc")

repName <- Vectorize(function(x){
  if(x != "Unique"){
    x = as.numeric(x)
    if(x < 10){
      x = paste0("0", as.character(x))
    }else{
      x = as.character(x)
    }
  }
  return(paste0("Rep" , x))
})
dataa$Family <- repName(dataa$Family)
###################
library(dplyr)
library(tidyr)
library(corrplot)
efm <- dataa %>%
  filter(Organism == "Enterococcus faecium") %>%
  select(Family, Strain) %>%
  table()
M <- cor(t(efm))
corrplot(M, method = "circle")         


library(pheatmap)

pheatmap(efm)
colSums(efm)

pp <- function(dat, species)
{
  a <- dat %>%
    filter(Organism == species) %>%
    select(Family, Strain) %>%
    table()
  a<- as.matrix(a)
  
  a <- a[,order(colnames(a))]
  
  library("pheatmap")
  library("RColorBrewer")
  col.pal <- brewer.pal(9,"Blues")
  pheatmap(t(a), 
           cluster_row = F,
           cluster_cols = T,
           color = col.pal, 
           fontsize = 7,
           fontsize_row=7, 
           fontsize_col = 1.6,
           gaps_col=50,
           cellheight = 8)
  
  pheatmap(t(a[,-dim(a)[2]]), 
           cluster_row = F,
           cluster_cols = T,
           color = col.pal, 
           fontsize = 7,
           fontsize_row=7, 
           fontsize_col = 1.6,
           gaps_col=50,
           cellheight = 8)
}

pp(dataa, "Enterococcus faecalis")
pp(dataa, "Enterococcus faecium")

Dup_ef_strain <- dataa %>% filter(Organism == "Enterococcus faecalis") %>%
  select(Strain, Assembly) %>%
  group_by(Strain) %>%
  distinct(Assembly) %>%
  summarise(n = n()) %>%
  filter(n != 1) %>% 
  pull(Strain)

dataa %>% filter(Organism == "Enterococcus faecalis") %>%
  filter(Strain %in% Dup_ef_strain) %>%
  select(Family, Strain, Assembly, Acc) %>%
  arrange(Strain)


dataa %>% filter(Organism == "Enterococcus faecalis") %>%
  filter(!is.na(Assembly)) %>%
  write.table(file = "EF_rep_info.tsv", sep = "\t", row.names = F)



Dup_efm_strain <- dataa %>% filter(Organism == "Enterococcus faecium") %>%
  select(Strain, Assembly) %>%
  group_by(Strain) %>%
  distinct(Assembly) %>%
  summarise(n = n()) %>%
  filter(n != 1) %>% 
  pull(Strain)

dataa %>% filter(Organism == "Enterococcus faecium") %>%
  filter(Strain %in% Dup_efm_strain) %>%
  select(Family, Strain, Assembly, Acc) %>%
  arrange(Strain)


dataa %>% filter(Organism == "Enterococcus faecium") %>%
  filter(!is.na(Assembly)) %>%
  write.table(file = "EFM_rep_info.tsv", sep = "\t", row.names = F)
