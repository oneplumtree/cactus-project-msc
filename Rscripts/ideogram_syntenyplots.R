
#install.packages("RIdeogram")
#Asynt found at https://github.com/simonhmartin/asynt

require(RIdeogram)
source("asynt.R")
library(dplyr)

### Whole Ideogram for Genome ----

#Make the karyotypes
cactus_karyotype <- read.table("QSchromsizes.txt", sep = "\t", header = T, stringsAsFactors = F)

sep1 <- strsplit (cactus_karyotype[,1],split = "_")
cactus_karyotype$Chr <- sapply (sep1,"[[",3)

renaming <- read.table("renamekaryotype.txt", header = T, stringsAsFactors = F)

cactus_karyotype$Chr <- renaming$new_karyotype[match(cactus_karyotype$Chr, renaming$old_karyotype)]

cactus_karyotype <- cactus_karyotype %>% arrange(Chr)

#Upload the gene densities - gff3 files generated from gmap, QS chromsizes generated from samtools faidx
gene_densitysmall <- GFFex(input = "singlelinesmallQSjour.gff3", karyotype = "QSchromsizes.txt", feature = "gene", window = 1000000)
gene_densitybig <- GFFex(input = "singlelinebigQSjour.gff3", karyotype = "QSchromsizes.txt", feature = "gene", window = 1000000)

gene_density <- rbind(gene_densitysmall, gene_densitybig)

sep1 <- strsplit (gene_density[,1],split = "_")
gene_density$Chr <- sapply (sep1,"[[",3)

gene_density$Chr <- renaming$new_karyotype[match(gene_density$Chr, renaming$old_karyotype)]

#alkaloid locations were found with blastn and gmap
alkaloids1 <- read.table("alkaloidideogramwholenew.txt", sep = "\t", header = T, stringsAsFactors = F)
alkaloids1 <- alkaloids1[, c(1:5,7)]

alkaloids1$Chr <- renaming$new_karyotype[match(alkaloids1$Chr, renaming$old_karyotype)]
alkaloids1$Chr <- as.character(alkaloids1$Chr)

#create the ideogram
ideogram(karyotype = cactus_karyotype, overlaid = gene_density, label = alkaloids1, label_type = "marker", output = "totalchromosome.svg", Lx = 20, Ly = 25 )
convertSVG("totalchromosome.svg", device = "png")

###Start of various synteny codes. 
#Note that the majority of the process is the same with the exception of manually fixing the karyotype files and renaming

###QS Small to Dragonfruit Synteny ----

alignments <- import.paf("mm20reorderedQSsmallpitaya.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSsmalldragonfruitkaryotypeshortened.txt", header = T)
condensedrename <- read.table ("renamepitaya.txt")
#subsetsynteny2$query <- gsub('\\D','', subsetsynteny2$query)
#subsetsynteny2$query <-gsub('0','', subsetsynteny2$query)
subsetsynteny2$query <- condensedrename$V2[match(subsetsynteny2[,4], condensedrename[,1])]
subsetsynteny2$reference <- condensedrename$V3[match(subsetsynteny2[,1], condensedrename[,2])]
#sep1 <- strsplit (subsetsynteny3[,1],split = "_")

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & (subsetsynteny3$query == 9 | subsetsynteny3$query == 8) ] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 9] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

#Loop to determine syntenic lengths and total lengths per chromosome set
lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}
#Took these numbers and copied them into excel for processing

sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])

chromo <- data.frame(subsetsynteny2[,c(1,4)])
table(chromo)
hist(table(chromo))

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")

###QS Dragonfruit to Small Inverse ----

alignments <- import.paf("mm20pitayareorderedQS.paf")

alignments2 <- subset(alignments, Rlen >= 200 & Qlen >= 200)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$query),]

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSsmalldragonfruitkaryotypeshortenedinverse.txt", header = T)

condensedrename <- read.table ("renamepitaya.txt")
subsetsynteny2$reference <- condensedrename$V2[match(subsetsynteny2[,1], condensedrename[,1])]

df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$query <- df$new_small[match(subsetsynteny2[,4], df[,4])]

subsetsynteny2$query <- condensedrename$V4[match(subsetsynteny2[,4], condensedrename[,2])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & (subsetsynteny3$query == 9 | subsetsynteny3$query == 8) ] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 9] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])



ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")


chromo <- data.frame(subsetsynteny3[,c(1,4)])
table(chromo)


###QS small and big ----


alignments <- import.paf("mm20genesRefsmallQSbig.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSkaryotypeshortened.txt", header = T)

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

#df <- data.frame(new = 1:11, old = 12:22)
df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$query <- df$new_big[match(subsetsynteny2[,4], df[,2])]
subsetsynteny2$reference <- df$new_small[match(subsetsynteny2[,1], df[,4])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)
#subsetsynteny3 <- subsetsynteny2[subsetsynteny2$query == 12,]


subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & (subsetsynteny3$query == 9 | subsetsynteny3$query == 2)] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

chromo <- data.frame(subsetsynteny3[,c(1,4)])

table(chromo)


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}



sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 2),3])



ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")

###QS big and small  ----


alignments <- import.paf("mm20ReflargeQSsmall.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSkaryotypeshortenedinverse.txt", header = T)

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

#df <- data.frame(new = 1:11, old = 12:22)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)
#subsetsynteny3 <- subsetsynteny2[subsetsynteny2$query == 12,]


subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]


df <- read.table("QSrenamebasic.txt", header = T)
subsetsynteny2$reference <- df$new_big[match(subsetsynteny2[,1], df[,2])]
subsetsynteny2$query <- df$new_small[match(subsetsynteny2[,4], df[,4])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

chromo <- data.frame(subsetsynteny2[,c(1,4)])
table(chromo)


lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 6),3])

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == (j)),3])
  print(savesum)
}


ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3)
convertSVG("chromosome.svg", device = "png")

#Dragonfruit to Large -----

alignments <- import.paf("mm20refDragonfruitBigQS.paf")

alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$query),]

sep1 <- strsplit (subsetsynteny2[,4],split = "_")
subsetsynteny2$query <- sapply (sep1,"[[",3)

dualkaryotype <- read.table("QSdragonfruitlargekaryotype.txt", header = T)

condensedrename <- read.table ("renamepitaya.txt")
subsetsynteny2$reference <- condensedrename$V2[match(subsetsynteny2[,1], condensedrename[,1])]

df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$query <- df$df_big[match(subsetsynteny2[,4], df[,2])]

#subsetsynteny2$query <- condensedrename$V4[match(subsetsynteny2[,4], condensedrename[,2])]


subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")


chromo <- data.frame(subsetsynteny3[,c(1,4)])
table(chromo)

lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)

for (j in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== 11 & lengthmatches[,2] == (j)),3])
  oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== 11 & !lengthmatches[,2] == (j)),3])
  
  print(savesum)
  print(oppositesum)
}


sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])


#Large Dragonfruit -----

alignments <- import.paf("mm20refQSlargedragonfruit.paf")


alignments2 <- subset(alignments, Rlen >= 500 & Qlen >= 500)
#alignments2 <- subset(alignments2, ref_data$seq_len >= 500000 & query_data$seq_len >= 500000)

synteny <- alignments2[,1:8]

synteny$fill <- "cccccc"
subsetsynteny2 <- synteny[grep("_scaffold", synteny$reference),]

sep1 <- strsplit (subsetsynteny2[,1],split = "_")
subsetsynteny2$reference <- sapply (sep1,"[[",3)

condensedrename <- read.table ("renamepitaya.txt")
subsetsynteny2$query <- condensedrename$V2[match(subsetsynteny2[,4], condensedrename[,1])]

df <- read.table("QSrenamebig.txt", header = T)
subsetsynteny2$reference <- df$df_big[match(subsetsynteny2[,1], df[,2])]

subsetsynteny2$reference <- as.numeric(subsetsynteny2[,1])
subsetsynteny2$query <- as.numeric(subsetsynteny2[,4])
str(subsetsynteny2)

subsetsynteny2[is.na(subsetsynteny2)] <- 0
subsetsynteny3 <- subset(subsetsynteny2, query > 0)

subsetsynteny3$fill <- "cccccc"

subsetsynteny3$fill[subsetsynteny3$reference == 9 & subsetsynteny3$query == 9] <- "f53c14"
subsetsynteny3$fill[subsetsynteny3$reference == 1 & subsetsynteny3$query == 1] <- "f58514"
subsetsynteny3$fill[subsetsynteny3$reference == 2 & subsetsynteny3$query == 2] <- "f5ec5a"
subsetsynteny3$fill[subsetsynteny3$reference == 3 & subsetsynteny3$query == 3] <- "a3f55a"
subsetsynteny3$fill[subsetsynteny3$reference == 4 & subsetsynteny3$query == 4] <- "27bf5e"
subsetsynteny3$fill[subsetsynteny3$reference == 5 & subsetsynteny3$query == 5] <- "2ae1d8"
subsetsynteny3$fill[subsetsynteny3$reference == 6 & subsetsynteny3$query == 6] <- "1e8ebf"
subsetsynteny3$fill[subsetsynteny3$reference == 7 & subsetsynteny3$query == 7] <- "9040cb"
subsetsynteny3$fill[subsetsynteny3$reference == 8 & subsetsynteny3$query == 8] <- "cb40b6"
subsetsynteny3$fill[subsetsynteny3$reference == 10 & subsetsynteny3$query == 10] <- "c23575"
subsetsynteny3$fill[subsetsynteny3$reference == 11 & subsetsynteny3$query == 11] <- "a71f33"


subsetsynteny3b <- subsetsynteny3[!subsetsynteny3$fill == "cccccc",]

chromo <- data.frame(subsetsynteny2[,c(1,4)])
table(chromo)

lengthmatches <- data.frame(subsetsynteny2[,c(1,4,8)])
sum(lengthmatches$matches)


for (j in 1:11){
print(j)
savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (j)),3])
oppositesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & !lengthmatches[,2] == (j)),3])

print(savesum)
print(oppositesum)
}


n=numeric(0)
list_mat <- list()
for (i in 1:11){
  n[i]=5^i
  m=numeric(0)
  m=matrix(data=0,nrow=n[i],ncol=n[i])
  
  
  for (j in 1:11){
    for (k in 1:11){
  print(j)
  savesum <- sum(lengthmatches[which(lengthmatches[,1]== (j) & lengthmatches[,2] == (k)),3])
  print(savesum)
    }
  }
  list_mat[[i]] <- savesum #Holding Matrix
}

sum(lengthmatches[which(lengthmatches[,1]== 8 & lengthmatches[,2] == 9),3])
sum(lengthmatches[which(lengthmatches[,1]== 9 & lengthmatches[,2] == 8),3])

dualkaryotype <- read.table("QSdragonfruitlargekaryotypeinverse.txt", header = T)

ideogram(karyotype = dualkaryotype, synteny = subsetsynteny3b)
convertSVG("chromosome.svg", device = "png")
