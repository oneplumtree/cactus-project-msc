rm(list=ls(all=TRUE))
library ("gplots") 
library("RColorBrewer")
library("ggplot2")

alkaloidTPM <- read.table("wholepeyotealkaloidTPM.txt", header = T)
alkaloidTPM2 <- as.matrix(alkaloidTPM[, -1])
rownames(alkaloidTPM2) <- alkaloidTPM[,1]


heatmap.2 (alkaloidTPM2, 
           col = brewer.pal(n = 6, name = "Reds"), 
           scale = "column", 
           margins=c(8,8), 
           trace = "none", 
           cexCol = 1.1,
           cexRow = 1.1,
           Colv = NA, 
           Rowv = NA,
           key = T
           #sepwidth=c(0.05,0.1),
           #sepcolor = "white", #line colour it separates the columns with
           # rowsep=1:nrow(alkaloidTPM2),
           # colsep=1:ncol(alkaloidTPM2)
)
#heatmap(alkaloidTPM2, scale="column", col = "Reds")

barplot(alkaloidTPM, main="Total Alkaloid TPM",
        xlab="Plant", ylab="TPM")

barplot(t(alkaloidTPM2), main="Total Alkaloid TPM",
        xlab="Genes", ylab="TPM", beside = TRUE, legend = T)


ggplot(data = alkaloidTPM2,
       mapping = aes(x = stuname, y = value)) + 
  geom_col(position = position_dodge())



###Stringtie TPM
library(tidyverse)
stringtieTPM <- read.table("stringtieTPM.txt", header = F, sep = "\t")

sep1 <- strsplit (stringtieTPM[,11],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5),split = "TPM")
stringtieTPM$tpm <- sapply (sep2,"[[",2)


widestringtie <- pivot_wider (stringtieTPM, names_from = V2, values_from = V3)
widestringtie <- as.data.frame(widestringtie)
rownames(widestringtie) <- widestringtie[,1]
widestringtie2 <- as.matrix(widestringtie[, -1])


heatmap.2 (widestringtie2, 
           col = brewer.pal(n = 6, name = "Reds"), 
           scale = "row", 
           margins=c(8,8), 
           trace = "none", 
           cexCol = 1.1,
           cexRow = 1.1,
           Colv = NA, 
           Rowv = NA,
           key = TRUE,
           legend = TRUE)
           #sepwidth=c(0.05,0.1),
           #sepcolor = "white", #line colour it separates the columns with
           # rowsep=1:nrow(alkaloidTPM2),
           # colsep=1:ncol(alkaloidTPM2)

#heatmap(alkaloidTPM2, scale="column", col = "Reds")


barplot(widestringtie2, main="Total Alkaloid TPM",
        xlab="Section", ylab="TPM")

par (mar= c(5,5,5,5))
barplot(widestringtie2, main="Stringtie TPM",
        xlab="Section", ylab="TPM", beside = TRUE, legend.text = TRUE,
        args.legend = list (x = "topright", bty = "n", inset = c(-0.15,0)))


par (mar= c(4,4,3,3))
barplot(t(widestringtie2), main="Stringtie TPM",
        xlab="Genes", ylab="TPM", beside = TRUE, legend = T, cex.names = 0.85)
args.legend(cex = 1.5)

ggplot(data = stringtieTPM,
       aes(x = V1, y = V3, fill = V2, group = V2)) + 
  geom_col(position = position_dodge())

p <- ggplot(data = stringtieTPM,
       aes(x = V2, y = V3, fill = V1, group = V1)) +
  xlab("Sections")+
  ylab("TPM")+
  geom_col(position = position_dodge()) +
  theme_classic()+
  scale_fill_grey()
p + labs (fill = "Genes")

##QS

library(tidyverse)
QSstringtieTPM <- read.table("QSstringtieTPM.txt", header = F, sep = "\t")

sep1 <- strsplit (QSstringtieTPM[,10],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",5),split = "TPM")
QSstringtieTPM$tpm <- sapply (sep2,"[[",2)

QSstringtieTPM$tpm <- as.numeric(QSstringtieTPM$tpm)
QSstringtieTPM <- as.data.frame(QSstringtieTPM[c(1,3,11)])
widestringtie <- pivot_wider (QSstringtieTPM, names_from = V1, values_from = tpm)
widestringtie <- as.data.frame(widestringtie)
rownames(widestringtie) <- widestringtie[,1]
widestringtie2 <- as.matrix(widestringtie[,-1])

QSstringtieTPM <- widestringtie2

heatmap.2 (widestringtie2, 
           col = brewer.pal(n = 6, name = "Reds"), 
           scale = "row", 
           margins=c(8,8), 
           trace = "none", 
           cexCol = 1.1,
           cexRow = 1.1,
           Colv = NA, 
           Rowv = NA,
           key = TRUE,
           legend = TRUE)
#sepwidth=c(0.05,0.1),
#sepcolor = "white", #line colour it separates the columns with
# rowsep=1:nrow(alkaloidTPM2),
# colsep=1:ncol(alkaloidTPM2)

#heatmap(alkaloidTPM2, scale="column", col = "Reds")


barplot(QSstringtieTPM, main="Total Alkaloid TPM",
        xlab="Section", ylab="TPM")

t()

par (mar= c(5,5,5,5))
barplot(QSstringtieTPM, main="Stringtie TPM",
        xlab="Section", ylab="TPM", beside = TRUE, legend.text = TRUE,
        args.legend = list (x = "topright", bty = "n", inset = c(-0.05,0)))


par (mar= c(4,4,3,6))
barplot(t(QSstringtieTPM), main="Stringtie TPM",
        xlab="Genes", ylab="TPM", beside = TRUE, legend.text = T, cex.names = 0.85,
        args.legend = list (x = "topright", bty = "n", inset = c(-0.15,0)))




check_overlap_length <- function (to_match,all_start,all_stop,query_start,query_stop){
  
  c1 <- (all_start > query_start & all_start < query_stop) #all_start falls within qstart + qstop
  c2 <- (all_start < query_start & all_start > query_stop) #all_start falls within qstart + qstop [reversed]
  c3 <- (all_stop > query_start & all_stop < query_stop)   #all_stop falls within qstart + qstop
  c4 <- (all_stop < query_start & all_stop > query_stop)   #all_stop falls within qstart + qstop [reversed]
  
  c5 <- (all_start < query_start & all_stop > query_start) #query start falls within all_start and all_stop
  c6 <- (all_start > query_start & all_stop < query_start) #query start falls within all_start and all_stop [reversed]
  c7 <- (all_start < query_stop & all_stop > query_stop)   #query stop falls within all_start and all_stop
  c8 <- (all_start > query_stop & all_stop < query_stop)   #query stop falls within all_start and all_stop [reversed]
  
  test1 <- to_match[(c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8),]
  
  sub_astart <- all_start[(c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8)]
  sub_astop <- all_stop[(c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8)]
  
  query_range <- query_start:query_stop
  
  
  if (nrow (test1) > 0){
    test1$overlap <- NA
    
    for (yyy in 1:nrow(test1)){
      vec1 <- sub_astart[yyy]:sub_astop[yyy]
      
      test1$overlap[yyy] <- sum (vec1 %in% query_range)
      
    }
    
  }
  
  test1
}

#install.packages("ape")
library(ape)


stringtiedata <- read.gff("Jour_wholejour.gtf", na.strings = c(".", "?"), GFF3 = FALSE)
transcriptsdata <- read.gff("gmapJourmappingsingleline.gff3", na.strings = c(".", "?"), GFF3 = TRUE)
transcriptsdata1 <- read.gff("williamsiijourgmapwhole.gff3", na.strings = c(".", "?"), GFF3 = TRUE)

sep1 <- strsplit (transcriptsdata[,9],split = ";")
sep2 <- strsplit (sapply (sep1,"[[",2),split = "Name=")
transcriptsdata$name <- sapply (sep2,"[[",2)

toMatch <- c("TRINITY_DN10959_c0_g1_i1",
             "TRINITY_DN10959_c0_g1_i4",
             "TRINITY_DN26_c0_g1_i5",
             "TRINITY_DN9001_c0_g1_i1",
             "TRINITY_DN15955_c0_g2_i1",
             "TRINITY_DN983_c0_g1_i3"
             )

sorted <- transcriptsdata[grep(paste(toMatch, collapse="|"),transcriptsdata$name),]

stringtie <- data.matrix(stringtiedata)
transcripts <- data.matrix(sorted)


sorted$results <- NA

for (i in 1:nrow(sorted)){
  
  sub_test <- stringtiedata[stringtiedata[i,1] == sorted[i,1],]	
  
  test1 <- check_overlap_length (sub_test,as.numeric (sub_test[,4]),as.numeric (sub_test[,5]),as.numeric (sorted[i,4]),as.numeric (sorted[i,5]))
  
  sorted$results[i] <- mean (as.numeric (test1$tpm))
  
  
}
#install.packages("tidyverse")
library(tidyverse)
sorted$stringtie_tpm <- NA
save$stringtie_tpm <- NA

for (i in 1:nrow(stringtiedata)){
  
  similar <- stringtiedata[which(stringtiedata[,1] %in% sorted[,1]),]
  
  test1 <- check_overlap_length (similar,as.numeric (similar[,4]),as.numeric (similar[,5]),as.numeric (stringtiedata[i,4]),as.numeric (stringtiedata[i,5]))
  
  save$stringtie_tpm[i] <- mean (as.numeric (test1$tpm))
  
  
}

stringtiedata[2,1]
transcriptsdata[2,1]

leveltranscripts <- transcriptsdata
transcriptsdata[,1] <- factor(transcriptsdata[,1], levels=levels(stringtiedata[,1]))

transcriptsdata3 <- leveltranscripts[!is.na(transcriptsdata[,1]),]


