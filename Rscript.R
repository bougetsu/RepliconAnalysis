setwd("~/Work/project/data/")

dat <- read.table("enter_geno_family_meta.tab", header = T, stringsAsFactors = F, sep = "\t",row.names = NULL, quote = "")

dim(dat)

head(dat)



res <- read.table("entero_genome.blastn@card_np_homolog_filtered.tab_uniq", header = T, stringsAsFactors = F, sep = "\t")

met <- read.table("All_geno_meta.tab", header = F, stringsAsFactors = F, sep = "\t")

colnames(met) <- c("Acc", "Name", "Org", "Date", "Strain", "mol", "isolate", "country", "plasmid")

str(met)

which(met$Strain == "")

norep <- met[-idx_have_rep,]


sort(unique(met$Strain))
sort(unique(dat$Strain))
idx <- unique(met$Strain) %in% unique(dat$Strain)
unique(met$Strain)[-idx]

met$Acc[idx_have_rep] %in% intersect(met$Acc, dat$Acc)


idx_have_rep <- met$Acc %in% dat$Acc
write.table(met[-idx_have_rep,], file = "NoRepGenmoe.tab", quote = F, row.names = F, sep = "\t")

q_rep <- dat[,1]

q_res <- res[,1]


both <- intersect(q_rep, q_res)

idx_rep <- q_rep %in% both
idx_res <- q_res %in% both

a1 <- dat[idx_rep, ]
a2 <- res[idx_res, c(1, 2, 6, 7)]

write.table(a1, file = "test_rep.tab", quote = F, row.names = F, sep = "\t")
write.table(a2, file = "test_res.tab", quote = F, row.names = F, sep = "\t")

aa <- rbind(a1[,c(1, 2, 3, 4)], a2[,c(1, 2, 6, 7)])


str(dat)

length(unique(dat$Strain))

length(unique(dat$Organism))

unique(dat$Acc)
table(dat$Organism)
a <- aggregate(Strain~Organism, dat, unique)
count <- apply(a, 1, function(x) length(x[2][[1]]))

a <- cbind(a, count)

sum(count)

a <- aggregate(Strain~Organism, a, length)

count <- a[,2][[1]]

dat1 <- dat[dat$Organism == "Enterococcus faecium",]
d1 <- data.frame(dat1$Strain, dat1$Family)

a<- as.matrix(table(d1))

a <- a[,c(1, 6, 7, 8, 9, 10,2, 3, 4, 5, 11)]



library("pheatmap")
library("RColorBrewer")
## if not installed, quickly add it as follows:
source("http://bioconductor.org/biocLite.R")
biocLite(c("RColorBrewer", "pheatmap"))

col.pal <- brewer.pal(9,"Blues")
pheatmap(t(a), 
         cluster_row = F,
         cluster_cols = T,
         color = col.pal, 
         fontsize = 7,
         fontsize_row=7, 
         fontsize_col = 7,
         gaps_col=50,
         cellheight = 20)

pheatmap(t(a[,-11]), 
         cluster_row = F,
         cluster_cols = T,
         color = col.pal, 
         fontsize = 7,
         fontsize_row=7, 
         fontsize_col = 7,
         gaps_col=50,
         cellheight = 20)


pheatmap(t(a), 
         cluster_row = F,
         cluster_cols = T,
         color = col.pal, 
         fontsize = 7,
         fontsize_row=7, 
         fontsize_col = 7,
         gaps_col=50,
         cellheight = 20,
         clustering_distance_cols = "euclidean")

pheatmap(t(a[,-11]), 
         cluster_row = F,
         cluster_cols = T,
         color = col.pal, 
         fontsize = 7,
         fontsize_row=7, 
         fontsize_col = 1.8,
         gaps_col=50,
         cellheight = 10)

pheatmap(t(a[,-8]), 
         cluster_row = F,
         cluster_cols = T,
         color = col.pal, 
         fontsize = 7,
         fontsize_row=7, 
         fontsize_col = 7,
         gaps_col=50,
         cellheight = 20,
         clustering_distance_cols = "correlation")

s <- a[,1]

head(dat)

length(unique(dat[which(dat$Organism == "Enterococcus faecium"),]$Strain))
length(unique(dat[which(dat$Organism == "Enterococcus faecalis"),]$Strain))
length(unique(dat[which(dat$Organism == "Enterococcus casseliflavus"),]$Strain))



pp <- function(dat, aa)
{
  dat1 <- dat[dat$Organism == aa,]
  d1 <- data.frame(dat1$Strain, dat1$Family)
  
  a<- as.matrix(table(d1))
  
  a <- a[,order(as.numeric(colnames(a)))]
  
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


pp(dat, "Enterococcus faecium")




alist <- read.csv("TableEnterococci.csv", header = T, stringsAsFactors = F)
str(alist)

dups <- alist$Strain[duplicated(alist$Strain)]

tdups <- alist[which(alist$Strain %in% dups),c(1, 3, 6)]

write.table(tdups, "DuplicatedStrainTable.csv", quote = F, row.names = F, sep = "\t")

acas <- read.table("ncbi-genomes-2018-03-19/acc_fasta2", header = T, sep = "\t", stringsAsFactors = F)

str(acas)
str(dat)
colnames(acas)[2] <- "Acc"
acas$Acc <- gsub("\\.[1-9]","",acas$Acc)
intersect(dat$Acc, acas$Acc)

dataa <- merge(dat, acas, by = "Acc")
datdup <- merge(dataa, tdups, by = "Assembly")


####check on Ass level0
drep <- data.frame(dataa$Assembly,dataa$Family)

colnames(drep) <- c("Assembly", "Family")
a <- table(drep)
a<- as.matrix(table(drep))

a <- a[,c(1, 8:13, 2:7, 14)]

trep <- cbind(Assembly = rownames(a), a, Sum = rowSums(a))
trep <- as.data.frame(trep)

##############
trep$Sum <- as.numeric(trep$Sum)
trep[which(trep$Sum == 18),]

#RepAss <- merge(trep, dataa[,c(5, 6, 8, 13)], by = "Assembly")
RepAss <- merge(trep, alist[,c(1, 3, 6)], by = "Assembly")



RepAssDup <- RepAss[which(RepAss$Strain %in% dups),]

RAD <- RepAssDup[order(RepAssDup$X.Organism.Name, RepAssDup$Strain),]
write.table(RAD, "AccRepStrain_duplated", quote = F, row.names = F, sep = "\t")
##########
largerhits <- RepAss[which(RepAss$Sum >= 12),]
largerhits <-largerhits[order(largerhits$Sum, largerhits$Strain),]
write.table(largerhits, "largerhits", quote = F, row.names = F, sep = "\t")
#########check on fasta level
drep <- data.frame(dataa$Acc,dataa$Family)

colnames(drep) <- c("Acc", "Family")
a <- table(drep)
a<- as.matrix(table(drep))

a <- a[,c(1, 8:13, 2:7, 14)]

trep <- cbind(Acc = rownames(a), a, Sum = rowSums(a))
trep <- as.data.frame(trep)

RepAss <- merge(trep, dataa[,c(1, 5, 6, 8, 13)], by = "Acc")
RepAssDup <- RepAss[which(RepAss$Strain %in% dups),]

RAD <- RepAssDup[order(RepAssDup$Organism, RepAssDup$Strain),]
idd <- which(duplicated(RAD$Acc))
RAD <- RAD[-idd,]

write.table(RAD, "NameRepStrain_duplated", quote = F, row.names = F, sep = "\t")

length(unique(dat$Strain))



##############remove dup gca
dupgca <- read.table("Dup_GCA")

##remove duplicated CGA and corresponding ACC FROM DATA
idx_gca <- which( !(alist$Assembly %in% dupgca[,1]) )

alist_nodup <- alist[idx_gca, ]

##find new dup
alist_nodup[which(alist_nodup$Strain  %in% alist_nodup[which(duplicated(alist_nodup$Strain)),3]), c(1, 3, 6)]

##merge with rep results

acas <- acas[which(acas$Assembly %in% alist_nodup$Assembly),]

dataa <- merge(dat, acas, by = "Acc")



####check on Ass level0
drep <- data.frame(dataa$Assembly,dataa$Family)

colnames(drep) <- c("Assembly", "Family")
a <- table(drep)
a<- as.matrix(table(drep))

a <- a[,c(1, 8:13, 2:7, 14)]

ssum <- rowSums(a)
trep <- cbind(Assembly = rownames(a), a)
trep <- as.data.frame(trep, stringsAsFactors = F)
trep <- data.frame(trep, Sum = ssum, stringsAsFactors = F)
##############
trep$Sum <- as.numeric(trep$Sum)
table(trep$Sum)

#RepAss <- merge(trep, dataa[,c(5, 6, 8, 13)], by = "Assembly")
RepAss <- merge(trep, alist_nodup[,c(1, 3, 6, 7, 8)], by = "Assembly")

RepAss$X.Organism.Name <- gsub("(\\w*) (\\w*) (.*)", "\\1 \\2",RepAss$X.Organism.Name)
RepAss$X.Organism.Name <- gsub("(\\w*) (\\w*)\\. (.*)", "\\1 \\2\\.",RepAss$X.Organism.Name)
RepAss <- RepAss[order(RepAss$X.Organism.Name, RepAss$Sum, decreasing = T),]

write.table(RepAss, "Assembly_Rep_nodup.xls", sep = "\t", quote = F, row.names = F)

largerhits <- RepAss[which(RepAss$Sum >= 10),]
largerhits <-largerhits[order(largerhits$Sum, largerhits$Strain),]
write.table(largerhits, "largerhits", quote = F, row.names = F, sep = "\t")




####get co-corrunce of res gene
res <- read.table("entero_genome.blastn@card_np_homolog_filtered.tab_uniq", header = T, stringsAsFactors = F, sep = "\t")
colnames(res)[1] <- "Acc"
##merge at Accession level, we want replicon and res on the same fragment

accrep <- data.frame(dataa$Acc,dataa$Family)

colnames(accrep) <- c("Acc", "Family")
a <- table(accrep)
a<- as.matrix(table(accrep))

a <- a[,c(1, 8:13, 2:7, 14)]

actrep <- cbind(Acc = rownames(a), a)
actrep <- as.data.frame(actrep, stringsAsFactors = F)
actrep <- data.frame(actrep, Sum = rowSums(a), stringsAsFactors = F)

actrep$Sum <- as.numeric(actrep$Sum)
table(actrep$Sum)

RepAcc <- merge(actrep, dataa[,c(1, 5, 6, 8, 13)], by = "Acc")

RepResAcc <- merge(RepAcc, res, by = "Acc")

res9 <- table(RepResAcc$X9, RepResAcc$subject_id)
res9 <- as.matrix(res9)[1:3,]

head(res9[,order(res9[3,], res9[2,], decreasing = T)])

RepResAcc$subject_id <- gsub(".*\\|(.*)$","\\1", RepResAcc$subject_id, perl = T)

write.table(RepResAcc, "RepResAccesionLevel.tab", quote = F, row.names = F, sep = "\t")

length(unique(RepResAcc$Strain))


RepRes <- read.table("RepResAccesionLevel.tab", header = T, sep = "\t", stringsAsFactors = F)

colnames(RepRes)

datr <- RepRes[, c(20, 21 )]

mat_res <- as.matrix(table(datr))
order(rownames(mat_res))
mat_res <- mat_res[order(rownames(mat_res)), ]
mat_res <- cbind(Assembly = rownames(mat_res), mat_res)
mat_res <- as.data.frame(mat_res, stringsAsFactors = F)

RepAss <- read.table("Assembly_Rep_nodup.xls", header = T, sep = "\t", stringsAsFactors = F )

mat_rep <- RepAss[, c(1:15, 17, 18)]

datrr <- merge(mat_rep, mat_res, by = "Assembly")

mat_rep$Assembly == rownames(mat_res)


a <- matrix(as.character(1:9), 3, 3)

mat_rep <- datrr[,2:15]
mat_res <- as.matrix(datrr[, 17:84])
class(mat_res) <- "numeric"

rr.cor <- cor(mat_rep, mat_res)

rr.cor <- rr.cor[-3, ]
heatmap(rr.cor)
library(pheatmap)
pheatmap(rr.cor)


library(entropy)

mi.plugin(mat_rep[,1], mat_res[,1], )

install.packages("infotheo")

library(infotheo)

mutinformation(mat_rep[,1], mat_res[,1])

dim(mat_rep)
dim(mat_res)

muinfo <- matrix(numeric(), 14, 68)

for(i in 1:14)
{
  for(j in 1:68)
  {
    muinfo[i, j] <- mutinformation(mat_rep[,i], mat_res[,j])
  }
}

rownames(muinfo) <- colnames(mat_rep)

colnames(muinfo) <- colnames(mat_res)

pheatmap(muinfo)


cor(c(1, 0, 0), c(0, 1, 0))

cor(c(1, 0, 0), c(0, 0, 0))
cor(c(1, 0, 0), c(0, 0, 1))
cor(c(1, 0, 0), c(1, 0, 0))

cor(c(1, 0, 0), c(2, 0, 0))

mutinformation(c(1, 0, 0), c(0, 1, 0))
mutinformation(c(1, 0, 0), c(1, 0, 0))

mutinformation(c(1, 0, 0), c(2, 0, 0))
dist(c(1, 0, 0), c(0, 1, 0))

############E.f E.fm
idx.efm <- which(datrr$X.Organism.Name == "Enterococcus faecium")
idx.ef <- which(datrr$X.Organism.Name == "Enterococcus faecalis")

rr.efm <- datrr[idx.efm,]


efm_rep <- rr.efm[,2:15]
efm_res <- as.matrix(rr.efm[, 18:85])
class(efm_res) <- "numeric"

efm.cor <- cor(efm_rep, efm_res)


efm.cor <- efm.cor[c(1, 2, 5, 8, 10, 12, 13, 14), ]

efm.cor <- efm.cor[, -which(is.na(colSums(efm.cor)))]
pheatmap(efm.cor)


############E.f E.fm
rr.ef <- datrr[idx.ef,]


ef_rep <- rr.ef[,2:15]
ef_res <- as.matrix(rr.ef[, 18:85])
class(ef_res) <- "numeric"

ef.cor <- cor(ef_rep, ef_res)


ef.cor <- ef.cor[-c(3, 8, 10, 12), ]

ef.cor <- ef.cor[, -which(is.na(colSums(ef.cor)))]
pheatmap(ef.cor)

############together
rr.eff <- datrr[c(idx.ef, idx.efm),]


eff_rep <- rr.eff[,2:15]
eff_res <- as.matrix(rr.eff[, 18:85])
class(eff_res) <- "numeric"

eff.cor <- cor(eff_rep, eff_res)


eff.cor <- eff.cor[-c(3), ]

eff.cor <- eff.cor[, -which(is.na(colSums(eff.cor)))]
pheatmap(eff.cor)

############pie chart statistic with Ef
RepAss <- read.table("Assembly_Rep_nodup.xls", header = T, sep = "\t", stringsAsFactors = F )
rep.efm <- RepAss[RepAss$X.Organism.Name == "Enterococcus faecium", ]
dim(rep.efm)
colSums(rep.efm[, c(2:15)])
colSums(rep.efm[, c(2:15)] != 0)


rep.ef <- RepAss[RepAss$X.Organism.Name == "Enterococcus faecalis", ]
dim(rep.ef)
colSums(rep.ef[, c(2:15)])
colSums(rep.ef[, c(2:15)] != 0)
#############################
###########
install.packages("devtools")
devtools::install_git("https://gitlab.com/sirarredondo/mlplasmids")
library(mlplasmids)
setwd("~/Work/project/data/ncbi-genomes-2018-03-19/")

files <- dir()[-c(1, 2)]

test1 <- files[1:20]
l <- list()
i <- 1
for(f in test1)
{
  l[[i]] <- plasmid_classification(path_input_file = f, species = 'Enterococcus faecium')
  i <- i+1
}
l
